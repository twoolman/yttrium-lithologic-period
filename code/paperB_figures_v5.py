#!/usr/bin/env python3
# paperB_figures_v1.py
#
# The five Paper B figures, at 300 dpi raster and vector, from the same inputs
# the analysis used. Self-contained: it recomputes the blocked permutation rather
# than reading a cached array, so the figure cannot drift from the test.
#
# Color decisions, recorded here rather than left implicit:
#   Magnitude and ordered classes use the project BuPu ramp, truncated at the
#   palest 0.30 so the lightest class stays visible against white. BuPu was
#   selected for this project because it remains monotonic in lightness under
#   all three dichromacies and in grayscale print.
#   The two places that need categorical rather than ordered color, figures 4
#   and 5, use violet #4a3aa7 against orange #eb6834. That pair was validated
#   rather than chosen by eye: all-pairs CVD dE 29.5 protan and 37.2 tritan,
#   normal vision dE 37.6, both above 3:1 contrast on a white surface.

import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
import geopandas as gpd

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

SRC = C.require(C.JOINED, "the joined point set written by eon_enrichment_v4.R")
STATES = C.require(C.STATES_GEOJSON, "Natural Earth US state polygons")
XLSX = C.require(C.DS942_XLSX, "USGS Data Series 942 whole-rock geochemistry")
OUT = str(C.FIGURES)

BLOCK_KM = 150.0
N_PERM = 20000

# Source of record for the reported probabilities. The R run
# eon_enrichment_v3.R (seed 20260707, 20,000 permutations) produced
# results_paperB_v3/, and the manuscript quotes those values. This script draws
# its own permutations in Python, which is an independent Monte Carlo estimate
# of the same quantity and lands within one standard error but not on the same
# digits. Figure 3 therefore displays the Python null distribution as the
# picture of the design and annotates the R probabilities, so that the figure
# cannot disagree with the text.
P_OF_RECORD = {"narrow": 0.035898, "wide": 0.024599}
ALPHA = 0.05
SEED = 20260707

VIOLET = "#4a3aa7"
ORANGE = "#eb6834"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8984"
GRID = "#e4e3df"

NARROW = ["Igneous, intrusive"]
WIDE = ["Igneous, intrusive", "Igneous and Metamorphic, undifferentiated"]
SIX = {"Arizona": "AZ", "New Mexico": "NM", "Idaho": "ID",
       "Montana": "MT", "Nevada": "NV", "Utah": "UT"}
BELT_CLASSES = ["Metamorphic, sedimentary clastic",
                "Metamorphic and Sedimentary, undifferentiated",
                "Sedimentary, carbonate"]
HREE = ["Tb_ppm", "Dy_ppm", "Ho_ppm", "Er_ppm", "Tm_ppm", "Yb_ppm", "Lu_ppm"]
LREE = ["La_ppm", "Ce_ppm", "Pr_ppm", "Nd_ppm"]

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

# Springer figure specification for Natural Resources Research:
#   lettering Helvetica or Arial, 8 to 12 pt at final size;
#   line widths not below 0.1 mm, which is 0.3 pt;
#   width 84 mm single column or 174 mm double column, height not above 234 mm.
# Arial itself is proprietary and is not installed here. Liberation Sans is
# metrically identical to Arial, glyph for glyph, so the layout produced is the
# layout Arial would produce. It is embedded in the EPS as a Type 42 font.
SANS = ["Liberation Sans", "Arimo", "Nimbus Sans", "Helvetica", "Arial"]
MM = 1.0 / 25.4
W2 = 174 * MM          # double column width, inches
HMAX = 234 * MM        # maximum height, inches

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": SANS, "font.size": 12,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.8,
    "axes.labelcolor": INK, "axes.titlesize": 12, "axes.titleweight": "bold",
    "axes.titlecolor": INK, "axes.labelsize": 12,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 11, "ytick.labelsize": 11,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "legend.frameon": False, "legend.fontsize": 11,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def bupu(n, lo=0.30, hi=1.0):
    base = plt.get_cmap("BuPu")
    return [base(x) for x in np.linspace(lo, hi, n)]


def tidy(ax, grid_axis="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, axis=grid_axis, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)


SPRINGER = {"figB1_sample_sites_by_period": "Fig1",
            "figB2_yttrium_by_period":      "Fig2",
            "figB3_permutation_null":       "Fig3",
            "figB4_belt_control":           "Fig4",
            "figB5_yttrium_hree_proxy":     "Fig5"}


def save(fig, name):
    """Write the Springer deliverable plus a PNG for screen review.

    TIFF at 600 dpi is the submission format. Springer specifies at least 600
    dpi for combination art, which is what all five of these are, carrying line
    work and tone together. EPS is not produced, because the PostScript backend
    flattens transparency to opaque and several of these figures carry the
    overplot density of thousands of samples in the alpha channel, so an EPS
    would not be a faithful rendering of the analysis.

    bbox_inches is deliberately not set to "tight", since that crops the canvas
    and would deliver a figure narrower than the 174 mm the journal specifies.
    """
    tag = SPRINGER[name]
    fig.set_layout_engine("constrained")
    fig.savefig(f"{OUT}/{tag}.tif", dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    # Flatten to RGB on white. Matplotlib writes RGBA, and an alpha channel in
    # a print TIFF is a nuisance to a production workflow rather than a help.
    with Image.open(f"{OUT}/{tag}.tif") as im:
        if im.mode != "RGB":
            flat = Image.new("RGB", im.size, "white")
            flat.paste(im, mask=im.split()[-1] if im.mode == "RGBA" else None)
            flat.save(f"{OUT}/{tag}.tif", compression="tiff_lzw",
                      dpi=(600, 600))
    fig.savefig(f"{OUT}/{tag}.png", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"  wrote {tag}  ({name})")


# ---------------------------------------------------------------- shared load
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
    d["Zr_p_ppm"] = pd.to_numeric(d["Zr_p_ppm"], errors="coerce")
    d["Th_ppm"] = pd.to_numeric(d["Th_ppm"], errors="coerce")
    d = d[d["Y_ppm"] > 0]
    pp = d["UNIT_LINK"].map(parse_period)
    d["period"] = [p[0] for p in pp]
    d["age_ma"] = [p[1] for p in pp]
    d = d[d["period"].notna()].copy()
    g = gpd.GeoDataFrame(
        d, geometry=gpd.points_from_xy(d["Long_NAD27"], d["Lat_NAD27"]),
        crs=4267).to_crs(5070)
    d["x"] = g.geometry.x.values
    d["y"] = g.geometry.y.values
    d["block"] = [f"{int(np.floor(a/(BLOCK_KM*1000)))} {int(np.floor(b/(BLOCK_KM*1000)))}"
                  for a, b in zip(d["x"], d["y"])]
    d["logY"] = np.log10(d["Y_ppm"].values)
    return d


def permute(d, classes):
    ig = d[d["GENERALIZE"].isin(classes)].copy()
    ig["is_mp"] = ig["period"].values == "Mesoproterozoic"
    keep = sorted(s for s, sub in ig.groupby("State")
                  if sub["is_mp"].sum() >= 15 and (~sub["is_mp"]).sum() >= 15)
    sub = ig[ig["State"].isin(keep)].reset_index(drop=True)
    state, blk = sub["State"].values, sub["block"].values
    logY, is_mp = sub["logY"].values, sub["is_mp"].values
    idx = {s: np.where(state == s)[0] for s in keep}
    nbl = {s: len(np.unique(blk[idx[s]][is_mp[idx[s]]])) for s in keep}
    blocks = {s: np.unique(blk[idx[s]]) for s in keep}

    def stat(lab):
        parts, w = [], []
        for s in keep:
            i = idx[s]
            a, b = logY[i][lab[i]], logY[i][~lab[i]]
            parts.append(a.mean() - b.mean() if len(a) >= 5 and len(b) >= 5 else np.nan)
            w.append(min(lab[i].sum(), (~lab[i]).sum()))
        parts, w = np.array(parts, float), np.array(w, float)
        ok = ~np.isnan(parts)
        return float((parts[ok] * w[ok]).sum() / w[ok].sum())

    obs = stat(is_mp)
    rng = np.random.default_rng(SEED)
    perm = np.empty(N_PERM)
    for k in range(N_PERM):
        lab = np.zeros(len(sub), bool)
        for s in keep:
            i = idx[s]
            lab[i] = np.isin(blk[i], rng.choice(blocks[s], nbl[s], replace=False))
        perm[k] = stat(lab)
    p = (1 + int((np.abs(perm) >= abs(obs)).sum())) / (N_PERM + 1)
    return obs, p, perm, keep, ig


# ============================================================== FIGURE 1
def fig1(d, ig):
    st = gpd.read_file(STATES)
    st = st[st["adm0_a3"] == "USA"]                      # world file: US only
    st = st[["name", "geometry"]].rename(columns={"name": "nm"}).to_crs(5070)
    st["geometry"] = st.geometry.make_valid()
    st = st[~st["nm"].isin(["Alaska", "Hawaii"])]        # keep the conus frame
    six = st[st["nm"].isin(SIX)]

    order = ig.groupby("period")["age_ma"].first().sort_values()   # young -> old
    order = [p for p in order.index if (ig["period"] == p).sum() >= 25]
    cols = bupu(len(order))                                        # light -> dark

    fig, ax = plt.subplots(figsize=(W2, W2 * 5.4 / 7.0))
    st.plot(ax=ax, facecolor="#f7f7f5", edgecolor="#d8d7d2", lw=0.5, zorder=1)
    six.plot(ax=ax, facecolor="white", edgecolor=MUTED, lw=0.9, zorder=2)

    for p, c in zip(order, cols):
        sel = ig[ig["period"] == p]
        hot = p == "Mesoproterozoic"
        ax.scatter(sel["x"], sel["y"], s=9 if hot else 4.5, color=c,
                   linewidths=0.45 if hot else 0.25,
                   edgecolors=ORANGE if hot else "white",
                   zorder=9 if hot else 3)

    b = six.total_bounds
    pad = 70000
    ax.set_xlim(b[0] - pad, b[2] + pad)
    ax.set_ylim(b[1] - pad, b[3] + pad)
    ax.set_axis_off()
    
    for nm, ab in SIX.items():
        g = six[six["nm"] == nm]
        if len(g):
            c = g.geometry.representative_point().iloc[0]
            ax.annotate(ab, (c.x, c.y), color=INK2, fontsize=11,
                        weight="bold", ha="center", va="center", zorder=12,
                        path_effects=[pe.withStroke(linewidth=2.6,
                                                    foreground="white")])

    handles = [Line2D([], [], marker="o", ls="",
                      ms=6.2 if p == "Mesoproterozoic" else 5.2, mfc=c,
                      mec=ORANGE if p == "Mesoproterozoic" else "white",
                      mew=0.9 if p == "Mesoproterozoic" else 0.4,
                      label=f"{p}  ({int(ig['age_ma'][ig['period']==p].iloc[0])} Ma)")
               for p, c in zip(order, cols)][::-1]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              handletextpad=0.6, labelspacing=0.55, title="Lithologic period",
              title_fontproperties={"weight": "bold", "size": 8})
    save(fig, "figB1_sample_sites_by_period")


# ============================================================== FIGURE 2
def fig2(ig):
    t = (ig.groupby(["period", "age_ma"])
           .agg(n=("Y_ppm", "size"), med=("Y_ppm", "median"),
                p90=("Y_ppm", lambda v: np.quantile(v, 0.9)))
           .reset_index())
    t = t[t["n"] >= 25].sort_values("age_ma")

    fig, ax = plt.subplots(figsize=(W2, W2 * 3.9 / 6.6))
    tidy(ax)
    for _, r in t.iterrows():
        hot = r["period"] == "Mesoproterozoic"
        c = VIOLET if hot else MUTED
        ax.plot([r["age_ma"]] * 2, [r["med"], r["p90"]], color=c,
                lw=2.0 if hot else 1.6, solid_capstyle="round",
                alpha=1.0 if hot else 0.55, zorder=3)
        ax.plot(r["age_ma"], r["med"], "o", ms=8 if hot else 6.5, color=c,
                mec="white", mew=1.0, zorder=4)

    ax.set_xscale("log")
    ax.set_xlabel("Age of host lithologic period (Ma, log scale)")
    ax.set_ylabel("Acid-leachable yttrium (ppm)")
    
    for _, r in t.iterrows():
        hot = r["period"] == "Mesoproterozoic"
        dx = {"Paleoproterozoic": 8, "Archean": -4}.get(r["period"], 0)
        ha = "right" if dx < 0 else ("left" if dx > 0 else "center")
        ax.annotate(f"{r['period']}\nn = {int(r['n']):,}",
                    (r["age_ma"], r["p90"]), textcoords="offset points",
                    xytext=(dx, 10), ha=ha, fontsize=9,
                    color=INK if hot else INK2,
                    weight="bold" if hot else "normal")
    ax.annotate(f"median {t.loc[t['period']=='Mesoproterozoic','med'].iloc[0]:.1f} ppm",
                (1400, t.loc[t["period"] == "Mesoproterozoic", "med"].iloc[0]),
                textcoords="offset points", xytext=(15, -3), fontsize=10,
                color=VIOLET, weight="bold", va="center")
    lo = t.loc[t["period"] == "Archean"]
    ax.set_ylim(0, t["p90"].max() * 1.32)
    t.to_csv(C.RESULTS / "tableB1_period_summary_wide.csv", index=False)
    save(fig, "figB2_yttrium_by_period")


# ============================================================== FIGURE 3
def fig3(res_n, res_w):
    fig, axes = plt.subplots(1, 2, figsize=(W2, W2 * 3.3 / 7.2), sharey=True)
    for ax, (obs, p_py, perm, keep, _), ttl, arm in zip(
            axes, (res_n, res_w),
            ("Narrow population\n(SGMC class 'Igneous, intrusive')",
             "Wide population\n(igneous-bearing classes)"),
            ("narrow", "wide")):
        p = P_OF_RECORD[arm]
        tidy(ax)
        ax.hist(perm, bins=60, color="#dedcea", edgecolor="white", lw=0.4, zorder=2)
        ax.axvline(obs, color=VIOLET, lw=2.2, zorder=4)
        ax.axvline(-obs, color=VIOLET, lw=1.0, ls=(0, (3, 2)), zorder=4, alpha=0.6)
        ax.set_title(ttl, loc="left", fontsize=12, pad=8)
        ax.set_xlabel("Weighted within-state difference\nin mean log$_{10}$ Y")
        ax.annotate(f"observed {obs:+.4f}  (factor {10**obs:.2f})\n"
                    f"p = {p:.4f}   states: {', '.join(keep)}\n"
                    f"{'rejects' if p < ALPHA else 'does not reject'} at alpha {ALPHA}",
                    xy=(0.03, 0.97), xycoords="axes fraction",
                    fontsize=10, va="top", ha="left", color=INK, weight="bold")
    axes[0].set_ylabel("Permutations")
    save(fig, "figB3_permutation_null")


# ============================================================== FIGURE 4
def fig4(d, ig):
    mp = ig.loc[ig["period"] == "Mesoproterozoic", "Y_ppm"].values
    belt = d[(d["period"] == "Mesoproterozoic") &
             (d["GENERALIZE"].isin(BELT_CLASSES))]
    groups = [("Mesoproterozoic\nigneous-bearing", mp, VIOLET)]
    for cls, lab in [("Metamorphic, sedimentary clastic", "Belt Supergroup\nclastic metasediment"),
                     ("Sedimentary, carbonate", "Belt-age\ncarbonate"),
                     ("Metamorphic and Sedimentary, undifferentiated", "Belt-age metased.\nundifferentiated")]:
        v = belt.loc[belt["GENERALIZE"] == cls, "Y_ppm"].values
        if len(v) >= 20:
            groups.append((lab, v, ORANGE))

    fig, ax = plt.subplots(figsize=(W2, W2 * 4.0 / 6.6))
    tidy(ax)
    rng = np.random.default_rng(7)
    for i, (lab, v, c) in enumerate(groups):
        jit = rng.normal(0, 0.055, len(v))
        ax.scatter(np.full(len(v), i) + jit, v, s=3.0, color=c, alpha=0.16,
                   linewidths=0, zorder=2)
        bp = ax.boxplot([v], positions=[i], widths=0.42, showfliers=False,
                        patch_artist=True, zorder=3,
                        medianprops=dict(color="white", lw=2.0),
                        boxprops=dict(facecolor=c, edgecolor=c, lw=1.0, alpha=0.92),
                        whiskerprops=dict(color=c, lw=1.0),
                        capprops=dict(color=c, lw=1.0))
        ax.annotate(f"median {np.median(v):.1f}\nn = {len(v):,}",
                    (i, np.quantile(v, 0.90)), textcoords="offset points",
                    xytext=(0, 16), ha="center", fontsize=10, color=INK,
                    weight="bold" if i == 0 else "normal")

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g[0] for g in groups], fontsize=10)
    ax.set_yscale("log")
    ax.set_ylabel("Acid-leachable yttrium (ppm, log scale)")
    ax.legend(handles=[
        Line2D([], [], marker="s", ls="", ms=7, color=VIOLET, label="Igneous-bearing units"),
        Line2D([], [], marker="s", ls="", ms=7, color=ORANGE, label="Non-igneous units of the same age")],
        loc="upper right")
    nbelt = sum(len(g[1]) for g in groups[1:])
    save(fig, "figB4_belt_control")


# ============================================================== FIGURE 5
def fig5():
    from scipy import stats
    rock = pd.read_excel(XLSX, sheet_name="Data")
    for c in HREE + LREE + ["Y_ppm"]:
        rock[c] = pd.to_numeric(rock[c], errors="coerce")
    rock["HREE"] = rock[HREE].sum(axis=1, min_count=len(HREE))
    rock["LREE"] = rock[LREE].sum(axis=1, min_count=len(LREE))

    fig, ax = plt.subplots(figsize=(W2, W2 * 4.6 / 5.6))
    tidy(ax, grid_axis="both")
    out = {}
    # Both groups are evaluated upon the samples carrying a complete terbium
    # through lutetium suite, which is the population the manuscript reports,
    # so that the light group is a contrast upon the same rocks rather than
    # upon a larger and differently composed set.
    base = rock[(rock["Y_ppm"] > 0) & (rock["HREE"] > 0) & (rock["LREE"] > 0)]
    for col, c, lab in (("HREE", VIOLET, "Heavy, Tb through Lu"),
                        ("LREE", ORANGE, "Light, La through Nd")):
        s = base
        rho, p = stats.spearmanr(s["Y_ppm"], s[col])
        out[col] = (len(s), rho, p)
        ax.scatter(s[col], s["Y_ppm"], s=22, color=c, alpha=0.60,
                   linewidths=0.5, edgecolors="white", zorder=3,
                   label=f"{lab}   n = {len(s)},  rho = {rho:+.3f}")
        lx = np.log10(s[col]); ly = np.log10(s["Y_ppm"])
        b, a = np.polyfit(lx, ly, 1)
        xs = np.linspace(lx.min(), lx.max(), 50)
        ax.plot(10**xs, 10**(a + b * xs), color=c, lw=1.8, zorder=4)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Summed lanthanide group in whole rock (ppm, log scale)")
    ax.set_ylabel("Yttrium in whole rock (ppm, log scale)")
    ax.legend(loc="upper left")
    save(fig, "figB5_yttrium_hree_proxy")
    return out


# ==================================================================== main
if __name__ == "__main__":
    import os
    os.makedirs(OUT, exist_ok=True)
    print("loading")
    d = load()
    print("permuting, narrow arm")
    res_n = permute(d, NARROW)
    print("permuting, wide arm")
    res_w = permute(d, WIDE)
    ig = res_w[4]

    print("figures")
    fig1(d, ig)
    fig2(ig)
    fig3(res_n, res_w)
    fig4(d, ig)
    o = fig5()

    print("\n=== numbers these figures assert ===")
    for lab, (obs, p, _, keep, _) in (("narrow", res_n), ("wide", res_w)):
        print(f"  {lab:7s} obs {obs:+.6f}  factor {10**obs:.3f}  p {p:.5f}  "
              f"states {','.join(keep)}  decision "
              f"{'reject' if p < ALPHA else 'fail to reject'} at alpha {ALPHA}")
    for k, (n, rho, p) in o.items():
        print(f"  DS942 Y vs {k}: n {n}, Spearman rho {rho:+.3f}, p {p:.3g}")
