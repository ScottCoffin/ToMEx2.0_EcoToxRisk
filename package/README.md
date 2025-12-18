# PSSDplusplus

R package for the ToMEx 2.0 probabilistic ecotoxicity workflow (Monte Carlo alignment + PSSD++).

## Installation

```r
# install.packages("devtools")
devtools::install_github("ScottCoffin/ToMEx2.0_EcoToxRisk", subdir = "package")
```

The workflow requires `ssdtools` version 0.3.7. The package includes a helper:

```r
PSSDplusplus::check_ssdtools_version()
```

## Data

The bundled ToMEx 2.0 ecotoxicity dataset is available as `PSSDplusplus::tomex2`
or via `data("tomex2", package = "PSSDplusplus")`.

## Quick start

```r
library(PSSDplusplus)

results <- run_pssd_reprex(
  sim = 10,
  n_sim = 10,
  workers = max(1, parallel::detectCores() - 1),
  overwrite_cache = FALSE
)

results$pnec_summary
```

## Walkthrough

Browse the end-to-end vignette that follows the workflow used in `PSSD_reprex.R`:

```r
vignette("pssdplusplus-walkthrough", package = "PSSDplusplus")
```

Outputs (PSSD cache and figures) default to temporary directories; supply `cache_dir` and `figure_dir` to persist them.
