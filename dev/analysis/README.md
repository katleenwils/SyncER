------------------------------------------------------------------------

# SyncER — manuscript analysis scripts

These scripts reproduce the benchmark and the cross-method comparison reported in the SyncER manuscript *(Wils & Ramisch, 2026*[^readme-1]*)*. They are **not** part of the package API and are **not installed** with the package — they live in the `dev/analysis/` folder of the source repository (excluded from the built package via `.Rbuildignore`) for transparency and reproducibility. To use them, clone the repository and run them from that folder, e.g.:

[^readme-1]: Wils, K, Ramisch, A. (2026) SyncER: Synchronicity testing for Event Records using tied age-depth models. *Scientific Reports.*

``` r
source("dev/analysis/method_comparison.R")
```

## Contents

| Script | What it does |
|-------------------|-----------------------------------------------------|
| `generate_dataset.R` | Builds the synthetic five-core benchmark dataset (`core_data/ (folder of CSVs)`) and a summary figure (`core_plots.pdf`) from prescribed age–depth models. Fully reproducible (fixed seeds). |
| `other_tests.R` | Reference implementations of the alternative synchronicity tests used for comparison: OxCal-style agreement/difference indices, Parnell et al. (2008) age differences, and the overlap coefficient. Functions only — sourced by `method_comparison.R`. |
| `method_comparison.R` | Runs the cross-method comparison across event deposits — SyncER's Synchronicity Score vs. the alternative tests — and writes the confusion-matrix / performance tables. |

## Running the comparison

`method_comparison.R` is self-contained: it takes **all input from the installed package** (`inst/extdata`) — the synthetic record data (`record_data_input`) and the Bacon age–depth model output for every core, both raw (`core*`) and synchronised (`core*_synced`).

``` r
source("dev/analysis/method_comparison.R")
```

Results (six CSVs — pairwise / combined / performance, each for raw and synchronised ages) are written to a fresh folder under `tempdir()`; set `results_dir` near the top of the script to keep them somewhere permanent. SyncER's synchronicity score is computed through the package itself (`compute_overall_synchronicity()`).

## Regenerating the dataset from scratch (optional)

``` r
source("dev/analysis/generate_dataset.R")
```

This writes `core_data/ (folder of CSVs)` and `core_plots.pdf` to `out_dir` (the working directory by default). To feed a freshly generated dataset through the full SyncER pipeline and produce your own Bacon output, follow `vignettes/workflow.Rmd`.

## Dependencies

Beyond SyncER itself, these scripts use: `overlapping`, `dplyr`, `tidyr` (comparison) and `ggplot2`, `gridExtra`, `cowplot`, `rintcal` (dataset generation). Install any that are missing before running.

## Reproducibility notes

- `generate_dataset.R` and the synchronicity score Monte Carlo use fixed seeds, so those steps are deterministic.
- Bacon age–depth modelling is **not** deterministic across runs, so the exact posterior ages depend on the specific Bacon run. The Bacon output used for the manuscript is bundled in `inst/extdata`; `method_comparison.R` reads it directly so its results are reproducible without re-running Bacon.
