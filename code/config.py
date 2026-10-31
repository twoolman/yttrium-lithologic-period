"""
Shared path configuration for the Paper B pipeline.

Every path in this project resolves from a single project root. The root is
taken from the HREE_ROOT environment variable when it is set, and otherwise
from the parent of the directory holding this file, which is the repository
root for a normal clone. Nothing in this project writes outside that root.

    data/      inputs, not tracked by git (see data/README.md)
    results/   outputs, not tracked by git
    results/figures/   the five manuscript figures

Set the root explicitly when the data live elsewhere:

    Linux/macOS   export HREE_ROOT=/path/to/project
    Windows       set HREE_ROOT=C:\\path\\to\\project
"""
import os
from pathlib import Path

ROOT = Path(os.environ.get("HREE_ROOT", Path(__file__).resolve().parent.parent))
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
GEOLOGY = DATA / "geology"

# Inputs, as named by their published releases.
NURE_CSV = DATA / "Reanalyzed_NURE-HSSRv8.csv"
DS942_XLSX = DATA / "1.4GaIntrsChemModData.xlsx"
DS942_GEOCHRON = DATA / "1.4GaIntrsGeochronData.xlsx"
STATES_GEOJSON = DATA / "states.geojson"

# Intermediate written by eon_enrichment_v4.R and consumed by the Python stages.
JOINED = RESULTS / "joined_geology_all6.csv.gz"

for _d in (DATA, RESULTS, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)


def require(path, what):
    """Fail loudly and usefully when an input is missing."""
    if not Path(path).exists():
        raise SystemExit(
            f"Missing input: {path}\n"
            f"  This is {what}.\n"
            f"  See data/README.md for where to obtain it, or set HREE_ROOT\n"
            f"  if your project directory is not {ROOT}."
        )
    return Path(path)
