# Lithologic period and yttrium enrichment in western US stream sediments

Analysis code for a study testing whether the lithologic period of intrusive
igneous terrain predicts acid-leachable yttrium in the stream sediment draining
it, across Arizona, Idaho, Montana, Nevada, New Mexico and Utah.

The headline result is that drainage over Mesoproterozoic intrusive terrain
carries a median acid-leachable yttrium concentration 2.03 times that of
drainage over intrusive terrain of all other periods, at an empirical two-sided
probability of 0.0246 under a blocked permutation test. The signal is magmatic
rather than chronological: a 2,072 sample Belt Supergroup control of comparable
age shows no elevation, and the Archean is the poorest period examined.

**Author:** Thomas A. Woolman
([ORCID 0000-0002-4120-2157](https://orcid.org/0000-0002-4120-2157))
**License:** MIT, see [LICENSE](LICENSE)

---

## What this repository contains

Five scripts, run in order. Each is heavily commented at the head with the
reasoning behind its design decisions, including the ones that did not work.

| Stage | Script | Language | Produces |
|---|---|---|---|
| 1 | `code/eon_enrichment_v4.R` | R | The spatial join, the period parse, the blocked permutation test, Tables B0 to B4, and the joined point set the later stages consume |
| 2 | `code/eon_perm_port_v2.py` | Python | An independent reimplementation of the permutation test, used to verify stage 1 rather than to extend it |
| 3 | `code/ds942_validate_v3.py` | Python | Validation of the yttrium proxy premise against USGS Data Series 942 whole-rock geochemistry |
| 4 | `code/ree_groundtruth_v2.py` | Python | Comparison of the predicted ground against the known rare earth occurrence inventory |
| 5 | `code/paperB_figures_v5.py` | Python | The five manuscript figures, at 600 dpi TIFF and 300 dpi PNG |

`code/config.R` and `code/config.py` hold every path the project uses. No other
file contains a path.

Stage 2 exists because stage 1 is the source of record and a result that
matters should not rest on one implementation. The Python port draws its own
Monte Carlo sample, so it lands within one standard error of the R
probabilities rather than on the same digits. It is a check, not a substitute.

## The reported numbers

The manuscript quotes `eon_enrichment_v4.R` at **seed 20260707** with **20,000
permutations** and a **150 km** block size. Those three constants are set at the
top of the script and are the only things that need to match for the reported
probabilities to reproduce exactly:

| Arm | Population | Observed ratio | Empirical two-sided p |
|---|---|---|---|
| Wide (primary) | Igneous, intrusive + Igneous and Metamorphic, undifferentiated | 2.033 | 0.024599 |
| Narrow (replication) | Igneous, intrusive | 1.651 | 0.035898 |

Alpha is 0.05. The script header explains at length why it was raised from
0.001, and why that revision is disclosed rather than quietly adopted.

## Quick start

```bash
git clone https://github.com/twoolman/yttrium-lithologic-period.git
cd yttrium-lithologic-period
# place the five public input datasets in data/ as described in data/README.md
Rscript code/eon_enrichment_v4.R      # stage 1, roughly 10 to 20 minutes
python3 code/eon_perm_port_v2.py      # stage 2
python3 code/ds942_validate_v3.py     # stage 3
python3 code/ree_groundtruth_v2.py    # stage 4
python3 code/paperB_figures_v5.py     # stage 5
```

Outputs land in `results/`, figures in `results/figures/`. Neither directory is
tracked by git.

[USER_GUIDE.md](USER_GUIDE.md) covers installation, the input data, what each
stage does, how to confirm the run reproduced, and what to change if you want
to apply the design to different ground.

## Data

No data are redistributed here. All five inputs are public government releases
and are fetched from their publishers, which keeps this repository small and
avoids relicensing material that is not mine to relicense.
[data/README.md](data/README.md) gives the source and the expected filename for
each one.

## Citing this

Please cite both the software and the paper. Machine-readable metadata is in
[CITATION.cff](CITATION.cff).

> Woolman, T. A. (2026). *Analysis code for lithologic period and yttrium
> enrichment in western United States stream sediments* (Version 1.0.0)
> [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

## Scope

This is the pipeline for one paper. It is research code written to answer a
specific question on a specific archive, not a general-purpose library. It is
documented and licensed so the result can be checked and the design reused, and
the reasoning is in the script headers because that is where it is useful.
