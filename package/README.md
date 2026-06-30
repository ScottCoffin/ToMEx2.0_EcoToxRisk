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

## Key functions

### Data alignment

| Function | Description |
|----------|-------------|
| `align_data()` | Align a toxicity dataset to food-dilution and tissue-translocation ERMs using a single set of parameters (deterministic). Core of the Monte Carlo workflow. |
| `MC_sim_align_parallel()` | Parallelised wrapper that runs `align_data()` across *n* Monte Carlo parameter draws. |
| `matrix_function()` | Generate a Latin Hypercube Sampling (Sobol') parameter matrix for Monte Carlo alignment. |

### Traditional (deterministic) SSDs

These functions fit model-averaged SSDs (log-logistic, log-normal, gamma, log-Gumbel) to pre-aligned data and return tier-specific hazard concentrations.

| Function | Description |
|----------|-------------|
| `make_tiered_SSDs(data, tier)` | Fit a tiered SSD. Returns a `tiered_SSD` S3 object with scalar HC estimates, plot data, and the raw ssdtools fit. |
| `print(result)` | Compact summary of species count and HCx estimate. |
| `ggplot2::autoplot(result, x_label, erm_label, base_size, point_size, label_size)` | Publication-ready SSD plot with model-averaged curve, CI ribbon, repelled species labels, and crosshair annotation at the actual HCx position. |
| `plotly_ssd(result, x_label, erm_label, lower_size, upper_size)` | Interactive Plotly SSD with hover tooltips and an arrow annotation pointing to the HCx on the curve. Requires the `plotly` package. |

**Tier behaviour** (tiers defined in [Mehinto et al. 2022](https://link.springer.com/article/10.1186/s43591-022-00033-3)):

| `tier` | Data filter | Collapse | HCx |
|--------|-------------|----------|-----|
| 1 | All endpoints | 25th %ile | HC5 lower confidence limit (most protective) |
| 2 | All endpoints | 25th %ile | HC5 estimate + 95% CI |
| 3 | Organism/Population, risk.13 ≠ 1 | Median | HC5 estimate + 95% CI |
| 4 | Organism/Population, risk.13 ≠ 1 | Median | HC10 estimate + 95% CI |

`make_tiered_SSDs()` returns a `tiered_SSD` S3 object. Call `print()` for a
compact summary and `ggplot2::autoplot()` to produce the SSD plot. The object
also exposes `predictions`, `collapsed_data`, `hc_data`, and `fit_dists` for
fully custom plots.

`make_tiered_SSDs()` expects the final dose column to be supplied as
`dose_new`; it does not apply assessment factors internally. In the example
below, `af.time` and `af.noec` are divided into the aligned exposure metric to
recreate the Mehinto et al. threshold preparation. If input doses are already
NOEC-equivalent, or if an unadjusted sensitivity run is intended, omit that
division before calling `make_tiered_SSDs()`.

Internal per-tier helpers (`SSD_function_t1`, `SSD_function_t2`, `SSD_function_t2_3`,
`SSD_function_t3_4`) remain in the package source for backward compatibility with
existing scripts but are no longer exported.

### Probabilistic SSDs (PSSD++)

| Function | Description |
|----------|-------------|
| `make_all_pSSDs()` | Fit PSSD++ models across all tier × environment × ERM combinations with parallelisation and caching. |
| `summarize_PNECs()` | Extract HC5/HC10 statistics from cached PSSD++ objects and pivot to a wide summary table. |
| `do.pSSD_mod()` | Per-species NOEC sampling engine with alignment-derived SDs (PSSD++). `UFt`/`UFdd` assessment-factor matrices are optional; use `apply_assessment_factors = FALSE` to treat them as 1. |
| `do.pSSD()` | Original per-species NOEC sampling engine (PSSD+, no alignment SD). `UFt`/`UFdd` are optional and default to unit multipliers when omitted. |

For PSSD++ workflows, `make_all_pSSDs()` and `run_pssd_reprex()` default to
`apply_assessment_factors = TRUE`, preserving the published correction using
`af.time` and `af.noec`. Set `apply_assessment_factors = FALSE` only for
sensitivity checks or data already expressed as NOEC-equivalent values.

### Utilities

| Function | Description |
|----------|-------------|
| `parameter_histograms_function()` | Visualise LHS parameter draws. |
| `run_pssd_reprex()` | End-to-end reproducible example pipeline. |
| `check_ssdtools_version()` | Verify that `ssdtools` 0.3.7 is installed. |

## Basic alignment → SSD workflow

```r
library(PSSDplusplus)
library(dplyr)

# 1. Filter ToMEx 2.0
tox <- tomex2 |>
  dplyr::filter(
    Group != "Bacterium", Group != "Plant",
    effect.metric != "HONEC",
    tier_zero_tech_f == "Red Criteria Passed",
    tier_zero_risk_f == "Red Criteria Passed",
    risk.13 != 0
  )

# 2. Align to ERMs using default parameters
aligned <- align_data(tox)

# 3. Prepare food-dilution surface-water subset
fd_water <- aligned |>
  dplyr::filter(
    environment %in% c("Freshwater", "Marine"),
    ingestible != "not ingestible",
    EC_env_v.particles.mL_ingest > 0,
    Group != "Algae",
    shape_f != "Not Reported", poly_f != "Not Reported"
  ) |>
  dplyr::mutate(
    dose_new = (EC_env_v.particles.mL_ingest / (af.time * af.noec)) * 1000
  ) |>
  tidyr::drop_na(dose_new) |>
  dplyr::filter(dose_new > 0)

# 4. Fit all four tiers — each returns a tiered_SSD S3 object
result_t1 <- make_tiered_SSDs(fd_water, tier = 1)  # HC5 LCL
result_t2 <- make_tiered_SSDs(fd_water, tier = 2)  # HC5 + CI
result_t3 <- make_tiered_SSDs(fd_water, tier = 3)  # HC5 + CI, org/pop
result_t4 <- make_tiered_SSDs(fd_water, tier = 4)  # HC10 + CI, org/pop

# Print compact summary
print(result_t3)

# Static publication plot (requires ggplot2 and ggrepel)
ggplot2::autoplot(result_t3, x_label = "particles/L", erm_label = "Food Dilution")

# Interactive plot with hover tooltips (requires plotly)
plotly_ssd(result_t3, x_label = "particles/L", erm_label = "Food Dilution")
```

For the full Monte Carlo uncertainty propagation workflow, see:

```r
vignette("pssdplusplus-walkthrough", package = "PSSDplusplus")
```

## Walkthrough

Browse the end-to-end vignette that follows the workflow used in `PSSD_reprex.R`:

```r
vignette("pssdplusplus-walkthrough", package = "PSSDplusplus")
```

Outputs (PSSD cache and figures) default to temporary directories; supply `cache_dir` and `figure_dir` to persist them.
