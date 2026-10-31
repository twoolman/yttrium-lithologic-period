## Shared path configuration for the Paper B pipeline.
##
## Every path in this project resolves from a single project root. The root is
## taken from the HREE_ROOT environment variable when it is set, and otherwise
## from the parent of the directory holding this file, which is the repository
## root for a normal clone. Nothing in this project writes outside that root.
##
##   data/      inputs, not tracked by git (see data/README.md)
##   results/   outputs, not tracked by git
##
## Set the root explicitly when the data live elsewhere:
##
##   Sys.setenv(HREE_ROOT = "C:/path/to/project")

.hree_this_file <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(f[1], mustWork = FALSE))
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE))
  }
  file.path(getwd(), "code", "config.R")
}

ROOT <- Sys.getenv("HREE_ROOT", unset = NA_character_)
if (is.na(ROOT) || !nzchar(ROOT)) {
  ROOT <- dirname(dirname(.hree_this_file()))
}
ROOT <- normalizePath(ROOT, winslash = "/", mustWork = FALSE)

DATA    <- file.path(ROOT, "data")
RESULTS <- file.path(ROOT, "results")
FIGURES <- file.path(RESULTS, "figures")
GEOLOGY <- file.path(DATA, "geology")

NURE_CSV   <- file.path(DATA, "Reanalyzed_NURE-HSSRv8.csv")
DS942_XLSX <- file.path(DATA, "1.4GaIntrsChemModData.xlsx")

## Intermediate written here and consumed by the Python stages.
JOINED <- file.path(RESULTS, "joined_geology_all6.csv.gz")

for (.d in c(DATA, RESULTS, FIGURES)) {
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)
}

hree_require <- function(path, what) {
  if (!file.exists(path)) {
    stop(sprintf(paste0("Missing input: %s\n  This is %s.\n",
                        "  See data/README.md for where to obtain it, or set\n",
                        "  HREE_ROOT if your project directory is not %s."),
                 path, what, ROOT), call. = FALSE)
  }
  invisible(path)
}

message(sprintf("[config] project root: %s", ROOT))
