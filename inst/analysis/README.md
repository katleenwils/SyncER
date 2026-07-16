---

editor_options: 
  markdown: 
    wrap: 72
---

# SyncER — manuscript analysis scripts

These scripts reproduce the benchmark and the cross-method comparison reported in the SyncER manuscript. They are **not** part of the package API — they are shipped inside the installed package for transparency and reproducibility. After installing SyncER you can locate this folder with:

``` r
system.file("analysis", package = "SyncER")
```

## Contents

| Script | What it does | Needs Bacon? |
|------------------|---------------------------|---------------------------|
| `generate_dataset.R` | Builds the synthetic five-core benchmark dataset (`core_data.xlsx`) and a summary figure (`core_plots.pdf`) from prescribed age–depth models. Fully reproducible (fixed seeds). | No |
| `other_tests.R` | Reference implementations of the alternative synchronicity tests used for comparison: OxCal-style agreement/difference indices, Parnell et al. (2008) age differences, and the overlap coefficient. Functions only — sourced by `method_comparison.R`. | No |
| `method_comparison.R` | Runs the cross-method comparison across event deposits — SyncER's Synchronicity Score vs. the alternative tests — and writes the confusion-matrix / performance tables. | No (uses bundled output) |

## Running the comparison

`method_comparison.R` is self-contained: it takes **all input from the installed package** (`inst/extdata`) — the synthetic record data (`record_data_input.xlsx`) and the Bacon age–depth model output for every core, both raw (`core*`) and synchronised (`core*_synced`). It does **not** re-run Bacon and does **not** run the SyncER synchronisation pipeline (that pipeline is demonstrated in the package vignette, `vignettes/workflow.Rmd`).

``` r
source(system.file("analysis", "method_comparison.R", package = "SyncER"))
```

Results (six CSVs — pairwise / combined / performance, each for raw and synchronised ages) are written to a fresh folder under `tempdir()`; set `results_dir` near the top of the script to keep them somewhere permanent. SyncER's synchronicity score is computed through the package itself (`compute_overall_synchronicity()`), so the comparison reflects the actual package method.

## Regenerating the dataset from scratch (optional)

``` r
source(system.file("analysis", "generate_dataset.R", package = "SyncER"))
```

This writes `core_data.xlsx` and `core_plots.pdf` to `out_dir` (the working directory by default). To feed a freshly generated dataset through the full SyncER pipeline and produce your own Bacon output, follow `vignettes/workflow.Rmd`.

## Dependencies

Beyond SyncER itself, these scripts use: `overlapping`, `dplyr`, `tidyr` (comparison) and `openxlsx`, `ggplot2`, `gridExtra`, `cowplot`, `rintcal` (dataset generation). Install any that are missing before running.

## Reproducibility notes

- `generate_dataset.R` and the synchronicity score Monte Carlo use fixed seeds, so those steps are deterministic.
- Bacon age–depth modelling is **not** deterministic across runs, so the exact posterior ages depend on the specific Bacon run. The Bacon output used for the manuscript is bundled in `inst/extdata`; `method_comparison.R` reads it directly so its results are reproducible without re-running Bacon.
