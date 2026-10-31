# Input data

None of these files are redistributed in this repository. All five are public
releases of the United States Geological Survey or of Natural Earth, and they
are better cited than copied: the publishers version them, and a stale copy
here would eventually disagree with the record without anyone noticing.

Download each into this directory under the filename shown, which is the name
the publisher ships. Nothing needs renaming.

Every stage checks its inputs before doing any work and fails with the missing
path and a description of what the file is, so a missing download is obvious
immediately rather than twenty minutes in.

---

## 1. Reanalyzed NURE-HSSR stream sediment geochemistry

**Filename:** `Reanalyzed_NURE-HSSRv8.csv` (about 25 MB)

The analytical archive the study rests on. Stream sediment samples from the
National Uranium Resource Evaluation Hydrogeochemical and Stream Sediment
Reconnaissance, reanalyzed by a near-total four-acid digestion.

Search the USGS ScienceBase catalog for "Reanalyzed National Uranium Resource
Evaluation Hydrogeochemical and Stream Sediment Reconnaissance". Version 8 is
what the published results used; a later version will shift the numbers
slightly.

> Smith, D. B., et al. (2018). *Geochemical and mineralogical data for soils
> and stream sediments of the conterminous United States.* U.S. Geological
> Survey.

The accompanying `NURE_Reanalysis_DataDictionary` CSV is worth downloading
alongside it. The pipeline does not read it, but it is what you will want when
a column name is unclear.

**Note on the digestion.** ALS ME-MS61L does not fully dissolve refractory
phases. Xenotime, which is yttrium phosphate and the mineral into which the
heavy lanthanides substitute, is among them. The measured yttrium therefore
understates the phase of interest, and the bias runs conservative.

---

## 2. State Geologic Map Compilation

**Filenames:** `geology/AZ.zip`, `geology/ID.zip`, `geology/MT.zip`,
`geology/NM.zip`, `geology/NV.zip`, `geology/UT.zip`
(50 to 300 MB each)

State-level geology polygons. The pipeline uses two fields, `UNIT_LINK` for the
map symbol and `GENERALIZE` for the lithologic class.

Download the six state packages from the USGS State Geologic Map Compilation
and place the archives in `geology/`. Stage 1 extracts each into
`geology/AZ/` and so on the first time it runs and reuses the extracted
directories afterward, so the archives can be deleted once extracted.

Each extracted directory must contain `AZ_geol_poly.shp` and its companion
files, under that naming.

> Horton, J. D. (2017). *The State Geologic Map Compilation (SGMC) geodatabase
> of the conterminous United States.* U.S. Geological Survey.

**Why the symbol is parsed rather than the age field read.** The published age
field is incomplete outside Arizona. The standard geologic map symbol carried
in `UNIT_LINK` is complete, so the period is recovered from it. That recovery
is itself one of the methodological contributions of the paper.

---

## 3. USGS Data Series 942 whole-rock geochemistry

**Filename:** `1.4GaIntrsChemModData.xlsx` (about 800 KB)

Whole-rock geochemistry of 1.4 Ga A-type granitoids, used to validate the
premise that yttrium tracks the heavy lanthanides rather than the light ones.
This is the independent test of the proxy, and it is the reason the proxy is
treated as an empirical claim rather than a mineralogical assumption.

The companion `1.4GaIntrsGeochronData.xlsx` is not read by the pipeline but is
part of the same release and is useful context.

> du Bray, E. A., Holm-Denoma, C. S., San Juan, C. A., Lund, K., Premo, W. R.,
> & DeWitt, E. (2015). *Geochemical, modal, and geochronologic data for 1.4 Ga
> A-type granitoid intrusions of the conterminous United States* (Data Series
> 942). U.S. Geological Survey.

---

## 4. State polygons

**Filename:** `states.geojson` (about 2 MB)

US state boundaries, used only to assign whole-rock samples to states and to
draw the base map in Figure 1.

Natural Earth admin-1 states and provinces, filtered to the United States, at
1:10m. Any equivalent state polygon layer works provided it carries a `name`,
`name_en` or `gn_name` field; the script takes whichever it finds.

---

## 5. Rare earth occurrence compilation

**Used by:** stage 4 only, `ree_groundtruth_v2.py`

The known rare earth occurrence inventory, used to compare the ground the model
flags against occurrences already on the record. Stage 4 runs its schema
reconnaissance and stops cleanly if the file is absent, so the first four
stages can be run without it.

> Bellora, J. D., Burger, M. H., Van Gosen, B. S., Long, K. R., Carroll, T. R.,
> Schmeda, G., & Giles, S. A. (2019). *Rare earth element occurrences in the
> United States* (Version 4.0). U.S. Geological Survey.

---

## Disk budget

The six geology archives dominate. Allow roughly 3 to 4 GB for the downloads
and the same again for the extracted shapefiles. Everything else totals under
30 MB.
