#!/usr/bin/env python3
# ds942_validate_v1.py
# Validation of Paper B premises against USGS Data Series 942
# (du Bray et al. 2015) whole-rock geochemistry of 1.4 Ga A-type granitoids.
#
# Three tests:
#   T1  Spatial assignment of the 1,162 rock samples to the six study states.
#   T2  Idaho Y/Zr decoupling: magmatic (present in rocks) or surficial (drainage only)?
#   T3  Rock-versus-drainage paired comparison over the same intrusions.

import numpy as np
import pandas as pd
import geopandas as gpd
from scipy import stats

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

XLSX = C.require(C.DS942_XLSX, "USGS Data Series 942 whole-rock geochemistry")
DRAIN = C.require(C.JOINED, "the joined point set written by eon_enrichment_v4.R")
STATES = C.require(C.STATES_GEOJSON, "Natural Earth US state polygons")
CRS_EQ = 5070
SIX = ["Arizona", "New Mexico", "Idaho", "Montana", "Nevada", "Utah"]
ABBR = {"Arizona": "AZ", "New Mexico": "NM", "Idaho": "ID",
        "Montana": "MT", "Nevada": "NV", "Utah": "UT"}
HREE = ["Tb_ppm", "Dy_ppm", "Ho_ppm", "Er_ppm", "Tm_ppm", "Yb_ppm", "Lu_ppm"]

def med(x):
    x = pd.to_numeric(pd.Series(x), errors="coerce").dropna()
    x = x[x > 0]
    return (np.median(x), len(x)) if len(x) else (np.nan, 0)

# ---------------------------------------------------------------- load rocks
rock = pd.read_excel(XLSX, sheet_name="Data")
rock = rock[np.isfinite(pd.to_numeric(rock["Latitude"], errors="coerce")) &
            np.isfinite(pd.to_numeric(rock["Longitude"], errors="coerce"))].copy()
rock_g = gpd.GeoDataFrame(
    rock,
    geometry=gpd.points_from_xy(rock["Longitude"], rock["Latitude"]),
    crs=4326).to_crs(CRS_EQ)

st = gpd.read_file(STATES)
namecol = [c for c in st.columns if c.lower() in ("name", "name_en", "gn_name")][0]
st = st[[namecol, "geometry"]].rename(columns={namecol: "StateName"})
st = st.to_crs(CRS_EQ)
st["geometry"] = st.geometry.make_valid()

rock_g = gpd.sjoin(rock_g, st, how="left", predicate="within").drop(
    columns=["index_right"], errors="ignore")
rock_g["St"] = rock_g["StateName"].map(ABBR)

print("=== T1  DS 942 rock samples by state ===")
print("  total georeferenced rock samples: %d" % len(rock_g))
tot = rock_g["StateName"].value_counts(dropna=False)
for k, v in tot.items():
    mark = "  *" if k in SIX else ""
    print("    %-22s %5d%s" % (str(k), v, mark))
insix = rock_g[rock_g["St"].notna()].copy()
print("  within the six study states: %d (%.1f%%)" %
      (len(insix), 100 * len(insix) / len(rock_g)))

# ------------------------------------------------- T2 Idaho Y/Zr in the rocks
print("\n=== T2  Y/Zr in rocks, by state ===")
insix["Y_num"] = pd.to_numeric(insix["Y_ppm"], errors="coerce")
insix["Zr_num"] = pd.to_numeric(insix["Zr_ppm"], errors="coerce")
insix["Th_num"] = pd.to_numeric(insix["Th_ppm"], errors="coerce")
insix["YZr"] = insix["Y_num"] / insix["Zr_num"]
insix["YTh"] = insix["Y_num"] / insix["Th_num"]

rows = []
for s in ["AZ", "NM", "ID", "MT", "NV", "UT"]:
    d = insix[insix["St"] == s]
    if not len(d):
        continue
    my, ny = med(d["Y_num"])
    mz, nz = med(d["Zr_num"])
    mr, nr = med(d["YZr"])
    mt, nt = med(d["YTh"])
    rows.append((s, len(d), my, mz, mr, nr, mt))
    print("  %-3s n=%4d   Y %6.1f   Zr %7.1f   Y/Zr %6.3f (n=%4d)   Y/Th %6.2f"
          % (s, len(d), my, mz, mr, nr, mt))

idr = insix.loc[insix["St"] == "ID", "YZr"].dropna()
oth = insix.loc[(insix["St"] != "ID") & insix["St"].notna(), "YZr"].dropna()
if len(idr) >= 5 and len(oth) >= 5:
    u, p = stats.mannwhitneyu(idr, oth, alternative="two-sided")
    print("\n  Idaho rock Y/Zr vs other five states (Mann-Whitney):")
    print("    ID n=%d median %.3f   others n=%d median %.3f   p = %.4g"
          % (len(idr), np.median(idr), len(oth), np.median(oth), p))
else:
    print("\n  Idaho rock Y/Zr: insufficient n (ID=%d) for a test" % len(idr))

# ------------------------------------- T3 rock vs drainage over same intrusions
print("\n=== T3  rock vs drainage, paired over the same intrusions ===")
dr = pd.read_csv(DRAIN, low_memory=False)
dr = dr.drop(columns=["index_right"], errors="ignore")
dr = dr[np.isfinite(pd.to_numeric(dr["Lat_NAD27"], errors="coerce")) &
        np.isfinite(pd.to_numeric(dr["Long_NAD27"], errors="coerce"))].copy()
dr_g = gpd.GeoDataFrame(
    dr, geometry=gpd.points_from_xy(dr["Long_NAD27"], dr["Lat_NAD27"]),
    crs=4269).to_crs(CRS_EQ)

RADIUS_M = 10000.0
near = gpd.sjoin_nearest(dr_g, insix[["Intrs_Name", "St", "geometry"]].rename(
    columns={"St": "rock_St"}), how="inner", max_distance=RADIUS_M,
    distance_col="dist_m")
near = near.sort_values("dist_m").drop_duplicates(subset=["Lab_ID"])
print("  drainage samples within %.0f km of a DS 942 rock site: %d"
      % (RADIUS_M / 1000, len(near)))

near["Y_num"] = pd.to_numeric(near["Y_ppm"], errors="coerce")
near["Zr_num"] = pd.to_numeric(near["Zr_p_ppm"], errors="coerce")
near["YZr"] = near["Y_num"] / near["Zr_num"]

pairs = []
for name, dsub in near.groupby("Intrs_Name"):
    rsub = insix[insix["Intrs_Name"] == name]
    if len(dsub) < 5 or len(rsub) < 5:
        continue
    ry, _ = med(rsub["Y_num"])
    dy, _ = med(dsub["Y_num"])
    rz, _ = med(rsub["YZr"])
    dz, _ = med(dsub["YZr"])
    pairs.append((name, len(rsub), len(dsub), ry, dy, dy / ry, rz, dz))

pf = pd.DataFrame(pairs, columns=["intrusion", "n_rock", "n_drain",
                                  "rock_Y", "drain_Y", "ratio",
                                  "rock_YZr", "drain_YZr"])
if len(pf):
    pf = pf.sort_values("ratio", ascending=False)
    print("  intrusions with >=5 rock and >=5 drainage samples: %d" % len(pf))
    print("  %-34s %5s %5s %8s %8s %7s %8s %8s"
          % ("intrusion", "nRk", "nDr", "rockY", "drainY", "ratio",
             "rkY/Zr", "drY/Zr"))
    for r in pf.itertuples(index=False):
        print("  %-34s %5d %5d %8.1f %8.1f %7.2f %8.3f %8.3f"
              % (r.intrusion[:34], r.n_rock, r.n_drain, r.rock_Y,
                 r.drain_Y, r.ratio, r.rock_YZr, r.drain_YZr))
    w, pw = stats.wilcoxon(np.log10(pf["drain_Y"]), np.log10(pf["rock_Y"]))
    print("\n  paired Wilcoxon on log10 Y, drainage vs rock:")
    print("    k=%d intrusions   median ratio %.3f   p = %.4g"
          % (len(pf), np.median(pf["ratio"]), pw))
    wz, pz = stats.wilcoxon(np.log10(pf["drain_YZr"]), np.log10(pf["rock_YZr"]))
    print("  paired Wilcoxon on log10 Y/Zr, drainage vs rock:")
    print("    k=%d   median drainage/rock Y/Zr ratio %.3f   p = %.4g"
          % (len(pf), np.median(pf["drain_YZr"] / pf["rock_YZr"]), pz))
    pf.to_csv(C.RESULTS / "ds942_rock_vs_drainage.csv", index=False)
else:
    print("  no intrusion met the paired sample-size threshold")

# ------------------------------------------- Y as HREE proxy, by state
print("\n=== Y as HREE proxy within the six states ===")
for c in HREE:
    insix[c] = pd.to_numeric(insix[c], errors="coerce")
insix["HREE_sum"] = insix[HREE].sum(axis=1, min_count=len(HREE))
sub = insix[insix["HREE_sum"].notna() & insix["Y_num"].notna() &
            (insix["HREE_sum"] > 0) & (insix["Y_num"] > 0)]
print("  complete-suite samples inside the six states: %d" % len(sub))
if len(sub) >= 8:
    rho, pr = stats.spearmanr(sub["Y_num"], sub["HREE_sum"])
    print("    Spearman rho(Y, sum HREE) = %+.3f   p = %.3g" % (rho, pr))
    print("    median Y / sum HREE = %.2f" % np.median(sub["Y_num"] / sub["HREE_sum"]))
for s in ["AZ", "NM", "ID", "MT", "NV", "UT"]:
    d = sub[sub["St"] == s]
    if len(d) >= 8:
        rho, pr = stats.spearmanr(d["Y_num"], d["HREE_sum"])
        print("    %-3s n=%3d  rho %+.3f  p %.3g  median Y/HREE %.2f"
              % (s, len(d), rho, pr, np.median(d["Y_num"] / d["HREE_sum"])))
    elif len(d):
        print("    %-3s n=%3d  (below threshold)" % (s, len(d)))

insix.drop(columns=["geometry"]).to_csv(
    C.RESULTS / "ds942_rocks_sixstate.csv", index=False)
print("\nwrote ds942_rocks_sixstate_v1.csv")
