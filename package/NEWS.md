PSSDplusplus 0.3.4
==================

CRAN readiness
--------------
- Fixed the bundled vignette, which previously ran the full Monte Carlo +
  PSSD++ pipeline against the entire ~13,000-row `tomex2` dataset and could
  hang for 10+ minutes during `R CMD check`. The vignette now runs
  end-to-end on the bundled `tomex2_mini` subset in about two minutes, with
  `n_sim` raised from 3 to 10 for a more representative demonstration.
- `MC_sim_align_parallel()`, `make_all_pSSDs()`, and `run_pssd_reprex()` no
  longer default to `parallel::detectCores() - 1` workers. Per CRAN policy,
  packages must not default to using multiple cores; these functions now
  default to sequential execution (`num_cores = 1` / `parallel = FALSE`,
  `workers = 1`). Parallel processing is still fully supported and strongly
  recommended for large jobs -- see the updated `@param` docs, examples, and
  README for how to opt in.
- Moved `sensobol`, `doSNOW`, `foreach`, `future`, `future.apply`,
  `progressr`, `tictoc`, `crayon`, `ggpubr`, `cols4all`, and `ggtext` from
  `Imports` to `Suggests`, cutting mandatory non-base install-time
  dependencies from 28 to 16. Removed the unused `doParallel` dependency
  entirely. Each
  formerly-mandatory package now degrades gracefully when not installed:
  - `doSNOW`/`foreach` (multi-core alignment) and `future`/`future.apply`/
    `progressr` (multi-core PSSD fitting/progress bars) fall back to
    sequential processing with a warning if requested but unavailable.
  - `crayon`/`tictoc` (colored console output/timing) silently degrade to
    plain messages.
  - `cols4all` (color palette) falls back to `scales::hue_pal()`.
  - `ggtext` (rich-text plot labels) falls back to plain `ggplot2` text/label
    geoms.
  - `ggpubr` (multi-panel arranged plots) falls back to returning the
    primary pSSD plot alone.
  - `sensobol` remains required specifically for `matrix_function()` (its
    Sobol'/LHS sampling engine) and now fails with an informative
    `install.packages()` message rather than an install-time hard
    dependency; `data.table` remains a hard `Import` because its `:=`
    assignment syntax requires NAMESPACE-level registration to work
    correctly and is not safely optional.
  Verified via the full test suite and a targeted smoke test that
  monkey-patches dependency detection to simulate each package being
  absent, for both the deterministic (`align_data`/`make_tiered_SSDs`) and
  full Monte Carlo + PSSD++ pipelines.
- Reworked the ToMEx2 upstream data-sync GitHub Action so it no longer
  auto-commits data, bumps the package version, or cuts a release on a
  schedule. It now performs a dry-run build/check against any upstream data
  change and emails the maintainer a status report and recommendation,
  leaving the decision to apply and release a data update to the maintainer.

PSSDplusplus 0.3.3
==================

Bug fixes
---------
- Fixed GitHub source-archive installs on Windows by excluding raw repository
  input files from generated archives. This keeps
  `remotes::install_github(..., subdir = "package")` focused on the package
  subtree and avoids extraction failures from deep/non-ASCII raw-data paths
  outside the package.

PSSDplusplus 0.3.2
==================

New features
------------
- Added `hc5_from_endpoints()`, an exported ML-inference entry point that
  accepts per-species linear-scale effect concentrations and computes an HCx
  distribution with the same internal PSSD++ sampler used by `make_all_pSSDs()`.

PSSDplusplus 0.3.1
==================

Minor update
------------
- Clarified package documentation for assessment factors (`af.time`,
  `af.noec`, `UFt`, and `UFdd`), including where they are applied in PSSD++
  sampling and where they are only shown in examples to recreate Mehinto et al.
  threshold preparation.
- Added `apply_assessment_factors` to `run_pssd_reprex()`, `make_all_pSSDs()`,
  `do.pSSD_mod()`, and `do.pSSD()`. The default remains `TRUE`; setting it to
  `FALSE` treats assessment factors as unit multipliers for sensitivity checks
  or data already expressed as NOEC-equivalent values.
- Made `UFt` and `UFdd` optional in the exported PSSD samplers. When omitted,
  the samplers use unit matrices so unadjusted NOEC sampling can run without
  assessment-factor inputs.

PSSDplusplus 0.2.0
==================

New features
------------
- Added `tomex2_mini`: a compact 10-species (5 freshwater, 5 marine) pre-filtered subset
  of the ToMEx 2.0 database for rapid testing and examples. Regeneration script is in
  `data-raw/create_tomex2_mini.R`.
- Added `PSSDplusplus_glossary` help page (`?PSSDplusplus_glossary`) with definitions for
  all acronyms and field-specific terms used in the package and manuscript (Coffin et al.
  2026, doi:10.1016/j.jhazmat.2025.141021).
- Vignette now uses `tomex2_mini` for a faster, self-contained end-to-end demonstration.

Documentation improvements
---------------------------
- Added `\dontrun{}` examples to all exported functions.
- Added `@description` and `@noRd` blocks to 20+ internal helper functions so error
  messages trace to meaningful source locations.
- Expanded `@param rmore_method` to explain the distributional assumption difference
  between `"step"` (trapezoidal, Wigger et al. 2020) and `"lognormal"` methods.
- Added `@details` note to `do.pSSD` / `do.pSSD_mod` clarifying that the dot follows
  legacy naming from Wigger et al. (2020), not R S3 dispatch conventions.
- Updated DESCRIPTION: improved abstract-style description with full method summary and
  DOI, corrected ORCID, removed placeholder contributor entry.
- Added glossary table to README and a full glossary section to the vignette.

Parallelism fixes
-----------------
- `MC_sim_align_parallel()` now automatically falls back to sequential processing when
  `n_sim <= 20`, where PSOCK cluster startup and data-serialisation overhead exceeds
  compute savings. `num_cores` is also capped to `n_sim`.
- `make_all_pSSDs()` internally caps `workers` to `nrow(combo_tbl)` (the number of
  tier × environment × ERM combinations); workers beyond that ceiling are never spawned.
- Vignette `workers` default changed from `detectCores() - 1` to
  `min(detectCores() - 1, 4)`, reflecting the 4-combo ceiling of the default workflow,
  with explanatory inline comments.
- Fixed duplicated parallel-setup logic across two vignette chunks; detection and derived
  variables (`workers`, `num_cores`, cache paths) are now computed exactly once.

PSSDplusplus 0.1.0
==================

Updates
-------
- Added sediment compartment support with kg-based dosing and alignment logic.
- Expanded Monte Carlo outputs to include sediment endpoints and conversions.
- Added LHS parameter toggles and sediment-aware parameter diagnostics.
- Improved robustness to sparse inputs (single-species handling, safe PSSD matrices).
- Added per-combination status summaries and clearer skip/error diagnostics.
- Updated PSSD and PNEC plotting labels to reflect particles/L vs particles/kg.

PSSDplusplus 0.0.2
==================

Updates
-------
- Added sediment alignment support using kg-based particle doses.
- Extended Monte Carlo outputs and thresholds to include sediment results.
- Added LHS parameter toggles and sediment-aware parameter plots.
- Fixed PSSD prep for single-species datasets and added sediment-aware PNEC labels.
- Fixed sediment alignment to respect exposure.route and derive particles/kg when only mass-based sediment doses are available.
- make_all_pSSDs now skips empty environment/ERM combinations cleanly to avoid downstream errors.
- Filled missing polydisperse length bounds when only a single length is reported to stabilize food-dilution alignment.
- PSSD generators now safely handle cases where all species lack data (returns empty matrices instead of error).
- make_all_pSSDs now reports a per-combination status summary and minimum species requirement.
- pSSD plot x-axis labels now reflect particles/L vs particles/kg based on environment.
