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
