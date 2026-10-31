# User guide

This guide covers installation, the input data, running the pipeline, checking
that it reproduced, and adapting the design to other ground. It assumes no
familiarity with the project beyond what is in the README.

---

## 1. Requirements

### R

R 4.2 or later, with three packages:

```r
install.packages(c("sf", "dplyr", "readr"))
```

`sf` is the only one that is ever awkward, because it binds to GDAL, GEOS and
PROJ. On Windows and macOS the CRAN binary carries them and nothing further is
needed. On Linux install the system libraries first:

```bash
# Debian / Ubuntu
sudo apt install libgdal-dev libgeos-dev libproj-dev libudunits2-dev
# Fedora / RHEL
sudo dnf install gdal-devel geos-devel proj-devel udunits2-devel
```

Stage 1 checks for all three packages before it does any work and stops with
the exact `install.packages()` call if any are missing.

### Python

Python 3.9 or later:

```bash
pip install numpy pandas geopandas scipy matplotlib openpyxl pillow
```

`geopandas` pulls in the same geospatial stack as `sf`. On a machine where that
proves difficult, conda handles it more gracefully than pip:

```bash
conda install -c conda-forge geopandas matplotlib scipy openpyxl pillow
```

### Fonts

Stage 5 sets the figure lettering to Arial or Helvetica, per Springer's
specification. Arial is proprietary and is frequently absent on Linux. The
script falls back through Liberation Sans and Arimo, both of which are
metrically identical to Arial glyph for glyph, so the layout produced is the
layout Arial would produce. Install either if the figures come out in the
matplotlib default:

```bash
sudo apt install fonts-liberation     # Debian / Ubuntu
sudo dnf install liberation-sans-fonts # Fedora / RHEL
```

### Hardware

Stage 1 holds roughly 32,000 points and six state geology polygon layers in
memory and runs 20,000 permutations twice. Budget 8 GB of RAM and 10 to 20
minutes. The Python stages are faster and lighter.

---

## 2. Where files live

Every path resolves from one project root. By default that is the repository
itself, so a normal clone needs no configuration.

```
yttrium-lithologic-period/
├── code/          the five stages plus config.R and config.py
├── data/          inputs you supply, not tracked by git
│   └── geology/   the six state geology archives
└── results/       everything the pipeline writes, not tracked by git
    └── figures/   Fig1 to Fig5 in TIFF and PNG
```

To keep the data somewhere else, set `HREE_ROOT` to a directory containing
`data/`:

```bash
export HREE_ROOT=/scratch/hree          # Linux, macOS
set HREE_ROOT=D:\hree                   # Windows cmd
$env:HREE_ROOT = "D:\hree"              # Windows PowerShell
```

```r
Sys.setenv(HREE_ROOT = "D:/hree")       # R, before sourcing
```

Both config files print the resolved root when they load, so the first line of
any run tells you where it is reading and writing.

**RStudio users:** open the repository as an RStudio project, or set the
working directory to the repository root, so that `code/config.R` can be found.
Stage 1 looks for the config in three places in order: beside the script when
run with `Rscript`, then `./config.R`, then `./code/config.R`. If none resolve
it stops with a message rather than guessing.

---

## 3. Input data

Five public releases. None is redistributed here. Download each into `data/`
under the filename given, which is the name the publisher ships.

See [data/README.md](data/README.md) for the download links and the citation
for each. In summary:

| File | Source | Approximate size |
|---|---|---|
| `Reanalyzed_NURE-HSSRv8.csv` | USGS reanalyzed NURE-HSSR geochemistry (Smith et al. 2018) | 25 MB |
| `geology/AZ.zip` … `geology/UT.zip` | USGS State Geologic Map Compilation (Horton 2017), six states | 50 to 300 MB each |
| `1.4GaIntrsChemModData.xlsx` | USGS Data Series 942 whole-rock geochemistry (du Bray et al. 2015) | 800 KB |
| `states.geojson` | Natural Earth admin-1 US state polygons | 2 MB |
| REE occurrence compilation | USGS rare earth occurrences (Bellora et al. 2019) | 5 MB |

Stage 1 unzips the geology archives into `data/geology/AZ/` and so on the first
time it runs, and reuses the extracted directories afterward. The archives can
be deleted once extracted.

Every stage checks its inputs before doing work and fails with the missing path,
what the file is, and a reminder about `HREE_ROOT`. Nothing fails silently or
half way through.

---

## 4. Running it

Run the stages in order. Stage 1 must complete before any Python stage, because
it writes `results/joined_geology_all6.csv.gz`, which all four consume.

```bash
Rscript code/eon_enrichment_v4.R
python3 code/eon_perm_port_v2.py
python3 code/ds942_validate_v3.py
python3 code/ree_groundtruth_v2.py
python3 code/paperB_figures_v5.py
```

Stages 2 through 5 are independent of one another and can run in any order.

If a stage fails part way, rerun that whole stage rather than resuming inside
it. The scripts are written to be rerun from the top and are not designed for
partial execution.

### Stage 1, in more detail

1. Reads the NURE-HSSR archive, keeping sediment samples with coordinates and a
   positive yttrium determination, discarding reference samples.
2. Projects to Albers equal area conic, EPSG 5070.
3. Joins each sample to the state geology polygon containing it, one state at a
   time, repairing invalid polygons where the published layer needs it.
4. Recovers the geological period from the map symbol in `UNIT_LINK`, matching
   longest symbols first so that `TR` is not read as `T` and `PP` and `IP` are
   not read as `P`. The published age field is incomplete outside Arizona,
   which is why the symbol is parsed rather than the field read.
5. Writes `results/joined_geology_all6.csv.gz` for the later stages.
6. Reports the Belt Supergroup control, which is declared and never tested
   against. Its purpose is to show that Mesoproterozoic age alone does not
   predict enrichment.
7. Runs both arms, wide and narrow, writing Tables B0 to B4 for each into
   `results/paperB/`.

The permutation permutes which **blocks** carry the period label, within state,
holding the per-state count of labelled blocks fixed. That preserves the
spatial structure and the provincial composition while destroying only the
association under test. Samples sit one to two kilometres apart, so a
conventional test on thousands of them would be enormously overpowered and
would reject almost any null.

---

## 5. Confirming the run reproduced

Open `results/paperB/permutation_test.txt`. With the shipped constants it
should report:

```
observed = +0.308203 (factor 2.033)   wide arm
permutations: 20000, two-sided empirical p = 0.02460
observed = +0.217641 (factor 1.651)   narrow arm
permutations: 20000, two-sided empirical p = 0.03590
```

`results/paperB/tableB0_primary_and_replication.csv` carries the same two rows
in machine-readable form, and the run log records the seed, the permutation
count and the alpha actually used.

Stage 2 is the independent check. It validates itself before it is trusted: it
must reproduce the observed statistics already on record, +0.217641 for the
narrow population and +0.308203 for the wide one. If it does not reproduce
both, it says so and nothing downstream of it should be used. Its probabilities
will differ from stage 1 in the third decimal, because it draws its own Monte
Carlo sample. Agreement within about one standard error is the pass condition,
not digit-for-digit equality.

Figure 3 displays the Python null distribution as the picture of the design and
annotates the R probabilities, so the figure cannot disagree with the text.

### If the numbers differ

- **Seed, block size or permutation count changed.** The three constants at the
  top of stage 1 are the ones that matter. Everything else is derived.
- **A geology layer was revised.** SGMC is versioned and the publisher has
  reissued state layers. A changed polygon boundary moves samples between
  units.
- **The NURE archive was reissued.** Version 8 is what the paper used.
- **`sf` or GEOS changed how invalid polygons are repaired.** This moves a small
  number of boundary samples. It should shift the fourth decimal, not the
  result.

The run log in `results/paperB/` records the versions in play, which is the
first place to look when a rerun disagrees with the paper.

---

## 6. Figures

Stage 5 writes `Fig1.tif` through `Fig5.tif` into `results/figures/` at 600 dpi,
flattened to RGB on white with the alpha channel stripped, LZW compressed, at
174 mm width against a 234 mm height ceiling. That is Springer's specification
for combination art. It also writes 300 dpi PNGs of the same figures for drafts
and slides.

To retarget another publisher, change `W2` and `HMAX` near the top of the
script. Everything downstream is expressed relative to those two constants, so
the figures re-lay themselves rather than scaling.

Line weights are set deliberately above the published minimum. Springer's floor
is 0.1 mm, which is 0.3 pt, and the strokes here run 0.45 to 0.8 pt, because
journal reduction to column width halves them again and a hairline that
survives on screen does not survive in print.

---

## 7. Adapting the design

The parts most likely to be useful elsewhere, in rough order of portability:

**The blocked permutation.** The function that permutes block labels within a
stratum is the reusable piece. It applies to any spatially autocorrelated
dataset where a categorical label is under test and a conventional test would
be overpowered by pseudoreplication.

**The map symbol parser.** `parse_period` in stage 1 recovers the geological
period from the standard geologic map symbol. It is a longest-match table and
extends to other symbol conventions by adding rows. This matters because
compiled national geology layers frequently carry an incomplete age field while
carrying a complete symbol field.

**The declared control.** The Belt Supergroup comparison is reported and never
tested against. Building a control into the design that can only embarrass the
hypothesis, and declaring it before the result, is the part of this design worth
copying regardless of subject.

**Changing the study region** means editing the `states` vector in stage 1 and
supplying the corresponding geology archives. The period parser and the
permutation are region-independent. The SGMC class names in the `arms` list are
specific to that compilation and would need remapping for a different source.

---

## 8. Reporting a problem

Open an issue with the stage that failed, the full console output including the
resolved project root printed by the config, your R or Python version, and your
`sf` or `geopandas` version. Most failures are a missing input or a geospatial
stack that did not install cleanly, and both are evident from those four things.
