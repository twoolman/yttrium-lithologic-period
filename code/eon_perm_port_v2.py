#!/usr/bin/env python3
# eon_perm_port_v1.py
#
# Verification port of the blocked permutation test in eon_enrichment_v1.R and
# eon_enrichment_v2.R, written to confirm that the stale alpha in
# results_paperB_v1/permutation_test.txt is the only thing wrong with that file.
#
# The port is validated before it is trusted: it must reproduce the observed
# statistics already on record, +0.217641 (factor 1.651) for the narrow
# population and +0.308203 (factor 2.033) for the wide one. If it does not
# reproduce both, nothing downstream of it is used.

import re
import numpy as np
import pandas as pd

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

SRC = C.require(C.JOINED, "the joined point set written by eon_enrichment_v4.R")
BLOCK_KM = 150.0
N_PERM = 20000
ALPHA = 0.05
SEED = 20260707

NARROW = ["Igneous, intrusive"]
WIDE = ["Igneous, intrusive", "Igneous and Metamorphic, undifferentiated"]

FINE = [("TR", "Triassic", 230.0), ("PP", "Pennsylvanian", 310.0),
        ("IP", "Pennsylvanian", 310.0), ("MZ", "Mesozoic undiff.", 170.0),
        ("PZ", "Paleozoic undiff.", 370.0), ("CZ", "Cenozoic undiff.", 30.0),
        ("Pg", "Paleogene", 45.0), ("Q", "Quaternary", 1.3),
        ("N", "Neogene", 12.0), ("T", "Tertiary", 30.0),
        ("K", "Cretaceous", 100.0), ("J", "Jurassic", 175.0),
        ("P", "Permian", 280.0), ("M", "Mississippian", 340.0),
        ("D", "Devonian", 400.0), ("S", "Silurian", 430.0),
        ("O", "Ordovician", 465.0), ("C", "Cambrian", 510.0),
        ("Z", "Neoproterozoic", 700.0), ("Y", "Mesoproterozoic", 1400.0),
        ("X", "Paleoproterozoic", 2050.0), ("W", "Archean", 3000.0),
        ("A", "Archean", 3000.0)]


def parse_period(code):
    if not isinstance(code, str) or len(code) < 3:
        return (None, np.nan)
    s = re.sub(r";.*$", "", code[2:])
    for sym, name, age in FINE:
        if s.startswith(sym):
            return (name, age)
    return (None, np.nan)


def load():
    d = pd.read_csv(SRC, low_memory=False)
    d = d[d["UNIT_LINK"].notna()].copy()
    d["Y_ppm"] = pd.to_numeric(d["Y_ppm"], errors="coerce")
    d = d[d["Y_ppm"] > 0]
    pp = d["UNIT_LINK"].map(parse_period)
    d["period"] = [p[0] for p in pp]
    d["age_ma"] = [p[1] for p in pp]
    d = d[d["period"].notna()].copy()
    import geopandas as gpd
    g = gpd.GeoDataFrame(
        d, geometry=gpd.points_from_xy(d["Long_NAD27"], d["Lat_NAD27"]),
        crs=4267).to_crs(5070)
    xy = np.c_[g.geometry.x.values, g.geometry.y.values]
    d["block"] = [f"{int(np.floor(a / (BLOCK_KM * 1000)))} "
                  f"{int(np.floor(b / (BLOCK_KM * 1000)))}" for a, b in xy]
    d["logY"] = np.log10(d["Y_ppm"].values)
    return d


def arm(d, classes, label):
    ig = d[d["GENERALIZE"].isin(classes)].copy()
    ig["is_mp"] = ig["period"].values == "Mesoproterozoic"

    keep = []
    for s, sub in ig.groupby("State"):
        if sub["is_mp"].sum() >= 15 and (~sub["is_mp"]).sum() >= 15:
            keep.append(s)
    keep = sorted(keep)

    sub = ig[ig["State"].isin(keep)].copy().reset_index(drop=True)
    state = sub["State"].values
    blk = sub["block"].values
    logY = sub["logY"].values
    is_mp = sub["is_mp"].values

    idx = {s: np.where(state == s)[0] for s in keep}
    n_mp_bl = {s: len(np.unique(blk[idx[s]][is_mp[idx[s]]])) for s in keep}
    blocks = {s: np.unique(blk[idx[s]]) for s in keep}

    def stat(lab):
        parts, w = [], []
        for s in keep:
            i = idx[s]
            a = logY[i][lab[i]]
            b = logY[i][~lab[i]]
            if len(a) < 5 or len(b) < 5:
                parts.append(np.nan)
            else:
                parts.append(a.mean() - b.mean())
            w.append(min(lab[i].sum(), (~lab[i]).sum()))
        parts = np.array(parts, float)
        w = np.array(w, float)
        ok = ~np.isnan(parts)
        return float((parts[ok] * w[ok]).sum() / w[ok].sum())

    obs = stat(is_mp)

    rng = np.random.default_rng(SEED)
    perm = np.empty(N_PERM)
    for k in range(N_PERM):
        lab = np.zeros(len(sub), bool)
        for s in keep:
            i = idx[s]
            pick = rng.choice(blocks[s], n_mp_bl[s], replace=False)
            lab[i] = np.isin(blk[i], pick)
        perm[k] = stat(lab)
    p = (1 + int(np.sum(np.abs(perm) >= abs(obs)))) / (N_PERM + 1)

    print(f"\n=== {label} ===")
    print(f"  states entering the test : {', '.join(keep)}")
    print(f"  n samples                : {len(sub)}")
    print(f"  observed                 : {obs:+.6f} (factor {10**obs:.3f})")
    print(f"  empirical two-sided p    : {p:.6f}  over {N_PERM} permutations")
    print(f"  decision at alpha {ALPHA:.3f}  : "
          f"{'REJECT the null' if p < ALPHA else 'fail to reject'}")
    print(f"  null mean {perm.mean():+.4f}  sd {perm.std(ddof=1):.4f}")
    for s in keep:
        i = idx[s]
        print(f"    {s}: n_mp {int(is_mp[i].sum()):4d}  mp_blocks {n_mp_bl[s]:3d}"
              f"  other_blocks {len(np.unique(blk[i][~is_mp[i]])):3d}"
              f"  ratio {np.median(10**logY[i][is_mp[i]]) / np.median(10**logY[i][~is_mp[i]]):.3f}")
    return obs, p, keep


if __name__ == "__main__":
    d = load()
    print("samples with a parsed period:", len(d))
    on = arm(d, NARROW, "ARM 1  narrow population, SGMC class 'Igneous, intrusive' only")
    ow = arm(d, WIDE, "ARM 2  wide population, igneous-bearing classes")

    print("\n=== PORT VALIDATION against the recorded R output ===")
    for got, want, lab in ((on[0], 0.217641, "narrow"), (ow[0], 0.308203, "wide")):
        ok = abs(got - want) < 1e-5
        print(f"  {lab:7s} observed: port {got:+.6f}  R {want:+.6f}  "
              f"{'MATCH' if ok else 'MISMATCH, do not use'}")
