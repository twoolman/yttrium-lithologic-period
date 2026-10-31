#!/usr/bin/env python3
# ree_groundtruth_v1.py
#
# PAPER B, RQ4. Ground-truth validation of the lithologic-period screen against
# known rare earth element occurrences.
#
# THE QUESTION THIS ANSWERS, AND WHY IT IS WORTH A SEPARATE RESEARCH QUESTION.
# Paper B as it stands establishes that a geochemical contrast exists between
# lithologic periods. It does not establish that the contrast finds anything.
# A screen that separates two populations of stream sediment is a statistical
# result; a screen that preferentially captures ground where rare earth
# occurrences are already known is an exploration tool. The second claim is
# harder to reject and is the one a reviewer at an economic geology journal
# will look for.
#
# THE FRAMEWORK, AND WHY THIS ONE.
# Mineral prospectivity mapping has a standard validation apparatus and this
# analysis uses it rather than inventing a bespoke test. The prediction-area
# (P-A) plot of Yousefi and Carranza (2015, Computers and Geosciences 79,
# 69-81) and the success-rate curve are what the field expects, they are what
# referees know how to read, and the editor-in-chief of the target journal is a
# co-author of the method. Using a home-made metric here would be a choice to
# be argued with rather than a choice to be checked.
#
# THE THREE CONFOUNDS THAT DECIDE WHETHER THE RESULT MEANS ANYTHING.
#
# (1) REPORTING BIAS. An occurrence is in the database because somebody looked
# and then wrote it down. People looked where they expected to find something,
# which is frequently the same reasoning the screen encodes. An apparent
# association can therefore be a record of past exploration judgment rather
# than evidence about the ground. This is controlled by restricting the whole
# analysis to the area with NURE sample coverage and by permuting within
# spatial blocks, so the null retains the clustering of both the samples and
# the occurrences.
#
# (2) TARGET MISMATCH. The screen predicts yttrium, which is a proxy for the
# HEAVY rare earths. A carbonatite-hosted light rare earth deposit of the
# Mountain Pass type is not what it is predicting, and counting such a deposit
# as a hit or a miss would both be wrong. Occurrences are therefore classified
# into HREE-relevant and LREE-dominant, and the primary test runs on the
# former, with the latter reported as a negative control: the screen SHOULD NOT
# capture LREE carbonatites above chance, and if it does, it is tracking
# something generic about mineralized ground rather than heavy rare earths.
#
# (3) SUPPORT. An occurrence is a point; a sediment sample integrates a
# catchment. Association is therefore tested at a stated distance rather than
# by point coincidence, and the distance is varied so the reader can see
# whether the result depends on the choice.
#
# STATUS: this script runs its schema reconnaissance and stops if the data file
# is absent, so it can be executed before the download to confirm the paths.

import os
import re
import sys
import glob
import zipfile
import numpy as np
import pandas as pd

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

PROJ = C.DATA
DRAIN = C.JOINED
OUT = C.RESULTS / "ree_groundtruth"

BLOCK_KM = 150.0
N_PERM = 20000
SEED = 20261008
RADII_KM = [5, 10, 20, 40]

SIX = ["AZ", "ID", "MT", "NM", "NV", "UT"]

# Minerals and deposit styles that carry the heavy rare earths plus yttrium.
# Xenotime is yttrium phosphate; the heavy lanthanides substitute into the
# yttrium site, which is the premise the whole paper rests on.
HREE_TERMS = [
    "xenotime", "yttrium", "yttr", "gadolinite", "fergusonite", "samarskite",
    "euxenite", "thalenite", "allanite", "zircon", "thortveitite",
    "ion adsorption", "ion-adsorption",
]
LREE_TERMS = [
    "bastnas", "bastnae", "monazite", "carbonatite", "cerium", "lanthan",
]


def find_ree_file():
    """Locate the USGS REE occurrence release wherever it landed."""
    pats = ["*REE*", "*ree_us*", "*rare*earth*", "*RareEarth*"]
    hits = []
    for p in pats:
        hits += glob.glob(os.path.join(PROJ, p))
        hits += glob.glob(os.path.join(PROJ, "**", p), recursive=True)
    hits = sorted(set(h for h in hits if os.path.isfile(h)))
    return hits


def load_occurrences(path):
    """Read whichever form of the release is present."""
    if path.lower().endswith(".zip"):
        zdir = os.path.join(OUT, "ree_unzip")
        os.makedirs(zdir, exist_ok=True)
        with zipfile.ZipFile(path) as z:
            z.extractall(zdir)
        inner = [f for f in glob.glob(os.path.join(zdir, "**", "*"), recursive=True)
                 if f.lower().endswith((".csv", ".shp", ".xlsx"))]
        if not inner:
            raise SystemExit(f"no csv, shp or xlsx inside {path}")
        path = sorted(inner)[0]
        print(f"  using {os.path.basename(path)} from the archive")
    if path.lower().endswith(".csv"):
        return pd.read_csv(path, low_memory=False, encoding="latin-1"), path
    if path.lower().endswith(".xlsx"):
        return pd.read_excel(path), path
    import geopandas as gpd
    g = gpd.read_file(path)
    return g, path


def report_schema(df):
    """The field names are not documented publicly, so discover them."""
    print("\n=== SCHEMA RECONNAISSANCE ===")
    print(f"  records: {len(df)}   columns: {len(df.columns)}")
    for c in df.columns:
        s = df[c]
        nn = s.notna().sum()
        samp = "; ".join(str(v)[:38] for v in s.dropna().unique()[:4])
        print(f"  {c:28s} n={nn:5d}  {samp[:96]}")
    # guess the fields the analysis needs
    guess = {}
    for role, pats in {
        "lat": [r"^lat", r"latitude", r"^y$", r"dec.*lat"],
        "lon": [r"^lon", r"^long", r"longitude", r"^x$", r"dec.*lon"],
        "state": [r"^state", r"st_?name", r"^adm"],
        "status": [r"status", r"develop", r"^type$", r"operation"],
        "commodity": [r"commod", r"ore", r"mineral", r"ree"],
        "deposit": [r"deposit", r"model", r"host", r"geolog"],
        "name": [r"^name", r"site", r"district", r"prospect"],
    }.items():
        for c in df.columns:
            if any(re.search(p, str(c), re.I) for p in pats):
                guess.setdefault(role, []).append(c)
    print("\n=== FIELD ROLE GUESSES, confirm before the analysis is trusted ===")
    for k, v in guess.items():
        print(f"  {k:10s} -> {', '.join(v)}")
    return guess


def classify(df, textcols):
    """HREE-relevant against LREE-dominant, from the free-text fields."""
    blob = df[textcols].astype(str).agg(" ".join, axis=1).str.lower()
    h = blob.apply(lambda s: any(t in s for t in HREE_TERMS))
    l = blob.apply(lambda s: any(t in s for t in LREE_TERMS))
    cls = np.where(h & ~l, "HREE", np.where(l & ~h, "LREE",
                   np.where(h & l, "BOTH", "UNCLASSIFIED")))
    return pd.Series(cls, index=df.index)


def main():
    os.makedirs(OUT, exist_ok=True)
    print("=== ree_groundtruth_v1.py ===")
    print("Looking for the USGS REE occurrence release in the project folder.")
    hits = find_ree_file()
    if not hits:
        print(f"""
  NOT FOUND in {PROJ}

  Needed: Bellora, J.D., Burger, M.H., Van Gosen, B.S., Long, K.R.,
  Carroll, T.R., Schmeda, G., and Giles, S.A., 2019, Rare Earth Element
  Occurrences in the United States (ver. 4.0, June 2019): U.S. Geological
  Survey data release, https://doi.org/10.5066/F7FN15D1

  Download USGS_REE_US_CSV.zip (about 148 KB) or the shapefile archive and
  drop it in the project root. This container cannot reach sciencebase.gov
  or mrdata.usgs.gov, the same restriction that applied to the SGMC
  shapefiles. Nothing else is required.

  Re-run this script once the file is present. It will print the schema
  first so the field roles can be confirmed before any test is computed.
""")
        return 1

    print(f"  candidates: {hits}")
    df, used = load_occurrences(hits[0])
    guess = report_schema(df)

    need = ["lat", "lon"]
    missing = [r for r in need if r not in guess]
    if missing:
        print(f"\n  cannot proceed: no field matched {missing}.")
        print("  Inspect the schema above and set the names explicitly.")
        return 1

    print(f"""
=== NEXT STAGE, held until the field roles above are confirmed ===

  1. Restrict occurrences to {', '.join(SIX)} and to the NURE coverage
     footprint, so the comparison is made only where the screen could have
     spoken.
  2. Classify each occurrence HREE-relevant or LREE-dominant from the mineral
     and deposit fields. The LREE set becomes the negative control.
  3. For radii of {RADII_KM} km, score every occurrence by the screen:
     (a) the binary Mesoproterozoic igneous-bearing indicator, and
     (b) the local acid-leachable yttrium percentile from the drainage data.
  4. Success-rate curve: cumulative share of HREE occurrences captured against
     cumulative share of area, with the 45 degree line as the random baseline.
  5. P-A plot after Yousefi and Carranza (2015): prediction rate against
     occupied area across thresholds, with the weight We taken at the
     intersection. We greater than zero marks an efficient criterion.
  6. Significance by the same blocked permutation the paper already uses, at
     {BLOCK_KM:.0f} km within state and {N_PERM} permutations, so the null
     preserves the clustering of occurrences as well as of samples.
  7. The negative control reported beside the primary result. If the screen
     captures LREE carbonatites as readily as HREE occurrences, it is finding
     mineralized ground generally and the heavy rare earth claim fails.
""")
    df.to_csv(os.path.join(OUT, "ree_occurrences_raw.csv"), index=False)
    print(f"  wrote {OUT}/ree_occurrences_raw.csv for inspection")
    return 0


if __name__ == "__main__":
    sys.exit(main())
