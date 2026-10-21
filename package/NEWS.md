PSSDplusplus 0.4.0
==================

This release fixes errors in how maximum ingestible particle size ("mouth
size") is assigned for the Food Dilution ERM, excludes in vitro data from every
species sensitivity distribution (SSD), and updates the bundled ToMEx 2.0 data
to database release v2.0.1 with corrected body lengths. **Food Dilution
alignments and HC5s change; Tissue Translocation results change where in vitro
records are removed or a corrected mouth size is below the translocation
limit.** Analyses run with earlier versions should be re-run. It supersedes the unreleased 0.3.5 / 0.3.6
development versions (see also the 0.3.5 notes below). Tracking:
ScottCoffin/ToMEx2.0_EcoToxRisk#6 and SCCWRP/ToMEx_AquaticOrganisms#5.

Bug fixes (errors in code)
--------------------------
- **No default 5 mm ingestible size.** `align_data()` set the maximum
  ingestible size of any record without a mouth size to the upper size limit
  of the assessment (`x2D_set`, 5000 um). For the volume-based Food Dilution
  ERM this inflated the mean ingestible particle volume by orders of
  magnitude and drove the aligned threshold (particles/L) down, often for
  the most sensitive species (e.g. *Moina macrocopa*, *Pseudechinus huttoni*,
  *Strombidium sulcatum*). Such records now get `x2M_ingest = NA`,
  `ingestible = "unknown"` and an `NA` Food Dilution ERM, so they are left out
  of Food Dilution SSDs; Tissue Translocation is unchanged. A message reports
  how many records and species are excluded. `missing_mouth = "x2D"` (also in
  `MC_sim_align_parallel()`) reproduces the old behaviour for comparison only.
  This affected every package workflow that aligns data: `align_data()`,
  `MC_sim_align_parallel()` and `run_pssd_reprex()`, and therefore any
  `do.pSSD()`, `make_tiered_SSDs()` or `hc5_from_endpoints()` result computed
  from their output.
- **Reported mouth sizes were overwritten.** The ToMEx 2.0 data preparation
  recomputed every mouth size from body length, replacing the mouth sizes
  reported in ToMEx 1.0 (Koelmans et al. 2020) for 13 species / life stages
  (2,769 records in v2.0.1; e.g. *Danio rerio* 1.53 -> 0.4 mm, *Daphnia
  magna* adult/juvenile 0.341 -> 0.115 mm, *Gammarus fossarum* 0.56 -> 0.125
  mm, *Ceriodaphnia dubia* early -> 0.036 mm; also *Hyalella azteca*,
  *Daphnia pulex*, *Artemia* spp. and *Hydra attenuata*). A reported mouth
  size (the record's, else one reported for the same species and life stage)
  is now always used; the body-length model is used only when none is
  reported.
- **Body-length model coefficients were ignored.** `align_data()` accepted
  `beta_log10_body_length_input` and `body_length_intercept_input` but never
  used them. It now predicts mouth size from `body.length.cm` (Jams et al.
  2020) for records without a reported size, using each simulation's sampled
  coefficients in `MC_sim_align_parallel()`. `tomex2$max.size.ingest.um` also
  holds the database's own predictions (`max.size.ingest.basis` "predicted
  (Jams et al. 2020)", 6,952 records); these are not treated as reported but
  re-predicted, so the sampled coefficients reach them (at the default
  coefficients the values are identical, so deterministic results do not
  change). Uncertainty is propagated only for predicted mouth sizes: reported
  (measured) mouth sizes are used as is in every simulation. The
  `mouth_size_source` argument of the 0.3.6 development version is kept so
  existing code runs, but `"body_length"` (which replaced reported mouth sizes
  with predictions) is deprecated: it warns and is treated as `"reported"`.
- **SSD confidence limits were not reproducible.** The manuscript pipeline
  fixes the random seed (`set.seed(99)`) before every ssdtools fit; the
  package did not, so bootstrapped limits, and therefore every Tier 1
  threshold (the lower confidence limit of the HC5), changed from run to run
  and did not match the published values (e.g. freshwater Tissue
  Translocation Tier 1: 274 vs 108 particles/L in one run). `make_tiered_SSDs()`
  and the internal SSD functions gain `seed = 99`; the fit runs under that
  seed and the caller's random number stream is restored afterwards
  (`seed = NULL` keeps the old behaviour). HC5 point estimates (Tiers 2-4)
  were not affected.
- **Deterministic SSDs depended on whether ssdtools was attached.** With
  ssdtools 0.3.7 (the version the package requires to reproduce the
  manuscript), fitdistrplus looks up each distribution's density function on
  the search path. Calling the package without `library(ssdtools)` silently
  dropped distributions from the model average (freshwater Food Dilution
  Tier 2 HC5: 10.8 instead of 59.6 particles/L). ssdtools is now attached
  before every fit. ssdtools 0.3.7 documents this requirement (its 1.0.0
  NEWS: "ssdtools had to be loaded to fit distributions"); the manuscript and
  the ToMEx 2.0 Shiny app attach it, so their results were not affected.
  ssdtools >= 2.0 is not affected either, but it averages distributions
  differently by default (the HCx of the model-averaged mixture distribution
  instead of the weighted mean of each distribution's HCx), which lowers
  Tier 2-4 HC5s by 0-47% for these data; with `est_method = "arithmetic"`
  it reproduces 0.3.7 to within 1%. Use 0.3.7 to reproduce the manuscript.
  PSSD++ does not use ssdtools.
- Parallel workers (`MC_sim_align_parallel(num_cores > 1)`,
  `make_all_pSSDs(parallel = TRUE)`) are new R sessions that did not inherit
  library paths set in the calling session, so they could load different
  versions of PSSDplusplus or its dependencies than the caller. They now use
  the caller's `.libPaths()`.
- `make_all_pSSDs()` (and `run_pssd_reprex()`) stopped with "Palette
  cols4all.area9d only supports up to 9 colors" when the data held more than 9
  taxonomic groups, as the full v2.0.1 data do in freshwater and marine
  waters. More than 9 groups now get evenly spaced hues; up to 9 keep the
  manuscript colours. Plot colours only.
- **In vitro data reached the SSDs.** No SSD code path filtered on in vitro /
  in vivo status, so in vitro records (rainbow trout anterior-kidney cell
  cultures, 10.1016/j.fsi.2021.04.014) entered the freshwater tissue
  translocation SSD (and, through the 5 mm default, the food dilution SSD).
  New exported `ssd_eligible()` keeps in vivo records only and is called by
  `process_environment_data()`, `run_pssd_reprex()`, and `make_tiered_SSDs()`
  and `hc5_from_endpoints()` (when a `vivo_f` column is present);
  `MC_sim_align_parallel()` output now carries `vivo_f`. The missing-mouth
  message no longer counts in vitro records. New exported `normalize_vivo()`
  maps every spelling of the field (the ToMEx 2.0 data preparation recoded
  lowercase "in vitro" to "In Vivo"; that spelling is not in the current data,
  so no record was affected).
- `matrix_function()` sampled `sim.body.length.intercept` with the slope's
  standard error (0.1376) instead of `se_body_length_intercept` (0.3222).
- `align_data()` set `x2M_ingest` / `x2M_trans` to `NA` for a record whose
  mouth size exactly equalled `x2D_set` / `upper.tissue.trans.size.um`
  (strict `<` / `>` comparisons); the boundary value is now used.
- Internal `process_environment_data()` returned surface-water Tissue
  Translocation thresholds in particles/mL (1000x too low) while Food Dilution
  was in particles/L; both are now particles/L. Only the internal Sobol
  wrapper used it; exported functions were not affected.
- `MC_sim_align_parallel()` dropped its `x1D_set` / `x2D_set` arguments
  (fixed in 0.3.5; see below).
- Data provenance: the ToMEx 2.0 data preparation overwrote the source of
  every ToMEx 1.0 body length (`body.size.source`, e.g. "fishbase.org",
  "study") with "reported". Sources are now kept. This did not change any
  value.

Data updates (refinements, not code errors)
-------------------------------------------
- `tomex2`, `tomex2_mini` and `inst/extdata/aoc_z_tomex2.RDS` are rebuilt from
  ToMEx 2.0 database release v2.0.1 (2026-06-30). Besides corrected body
  lengths (e.g. *Oryzias latipes*, Chisada et al. 2021: 3.5 -> 2.84 cm;
  *Strombidium sulcatum* now has a body length), v2.0.1 corrected a few
  records (e.g. *Hyalella azteca*, Au et al. 2015, relabelled early ->
  juvenile; exposure durations of some cladoceran NOECs 2 -> 4 d; one
  *Daphnia magna* LOEC removed; effect metrics of some *Moina macrocopa* and
  *Hyalella azteca* records). `tomex2` has 12,797 records.
- Body lengths of tested animals corrected at record level where the database
  used a length that was not that of the tested animals (evidence tier in
  `body.length.tier`):
  - *Pinctada margaritifera* (Gardon et al. 2018): 9.8 -> 5.9 cm (measured
    shell height; A).
  - *Mytilus galloprovincialis*: Capolupo et al. 2021, 15 -> 5.0 cm (A);
    Abidli et al. 2021, 15 -> 6.54 cm (A); Masia et al. 2021, 15 -> 6.54 cm,
    imputed (C; not reported in the study).
  - *Crepidula onyx* (Lo & Chan 2018, larvae): 6 -> 0.0624 cm (A).
  - *Tripneustes gratilla* (Kaposi et al. 2014, larvae): 10 -> 0.0478 cm
    (B; interpolated from another culture).
  - *Hyalella azteca* (Au et al. 2015): 0.78 -> 0.340 cm (10-d test) and
    0.564 cm (42-d test) (B). The reported mouth size (0.112 mm) still
    applies, so the alignment is unchanged.
  - *Moina macrocopa* (Kim & Rhee 2021): no length -> 0.1665 cm, the largest
    control mean read from the study's vector Fig. 2B (1562 / 1665 / 1591 um
    for F0 / F1 / F2; A). Previously excluded from food dilution (no length).
  - *Pseudechinus huttoni* (Richardson et al. 2021; MSc thesis 2020): no
    length -> 0.04777 cm, day-10 filtered-seawater control body length (A).
    The database coded this study In Vitro, but whole live larvae were
    exposed; its 21 records are corrected to In Vivo.
  - *Centropristis striata* juveniles (Stienbarger et al. 2021): no length ->
    3.95 cm, derived from 0.75 g with the FishBase length-weight relationship
    (C).
  - The DOI of Kaposi et al. 2014 is corrected (was "0.1021/es404295e").
  - *Hydra attenuata* (Murphy & Quinn 2018) checked: already freshwater.

  The corrections are listed, with sources, in `body_length_measured.csv` in
  the ToMEx 2.0 database repository (copy in this repository under
  `scripts/monte carlo/ref data/`).
- Body length is resolved per record as: row-level correction, else the
  record's own value, else a species x life-stage reference; it is never set
  to a default. Early or juvenile records that carry an adult size reference
  (e.g. *Crassostrea gigas* and *Pinctada margaritifera* larvae, *Amphibalanus
  amphitrite*, *Eriocheir sinensis* juveniles) keep that value but are flagged
  `body.length.basis = "adult reference (stage mismatch)"`.
- Records with no resolvable size are listed by the database build
  (`unresolved_body_length.csv`) instead of being given a default; in vitro
  records are listed in `excluded_in_vitro.csv`.

New features
------------
- **Checkpointing for long runs.** `MC_sim_align_parallel()` gains
  `checkpoint_dir` and `checkpoint_every`: simulations run in chunks, each
  finished chunk is saved atomically, and a re-run with the same inputs
  reloads finished chunks, so a crashed or interrupted run resumes where it
  stopped. `make_all_pSSDs()` already saved each PSSD++ combination as it
  finished; a cached combination is now reused only if it was computed from
  the same inputs (data, `sim`, `cv_uf`, `rmore_method`, `quantile_type`,
  `apply_assessment_factors`). Previously any file with the combination's
  name was reused, so a run with changed data or settings could silently
  return old results. Pass a persistent `base_cache_dir` for long runs; the
  default temporary directory does not survive a crash.
- New `tomex2` columns: `body.length.basis`, `body.length.tier`,
  `body.length.cm.database` (value before correction), `max.size.ingest.basis`
  ("reported (record)", "reported (species, life stage)", "override",
  "predicted (Jams et al. 2020)", "unresolved"), `jams_domain`,
  `vivo_f.database` (in vitro / in vivo coding before correction) and
  `ssd_exclusion_reason` ("in vitro", "vivo status missing", or "food dilution
  only: unresolved body length or mouth size").
- New exported `jams_domain()`: classifies body length against the
  calibration of the Jams et al. (2020) model ("within" >= 30 mm,
  "less_robust" 9 to < 30 mm, "outside_calibration" < 9 mm).
- `align_data()` returns `max.size.ingest.basis.used` and `jams_domain`, uses
  an optional `max.size.ingest.override.um` column (a measured gape or a
  documented ingestion ceiling) before the body-length model, and warns
  (`verbose = TRUE`) when mouth sizes are predicted for animals < 30 mm.
  `MC_sim_align_parallel()` reports the same count once per run.

Documentation and tests
-----------------------
- `?tomex2` documents the body-length, ingestible-size and in vivo fields
  and their provenance; `?align_data` documents the precedence of mouth sizes;
  `?ssd_eligible`, `?normalize_vivo` and `?jams_domain` are new; the glossary
  (`?PSSDplusplus-glossary`) describes reported-first mouth sizes, no default
  size, and in vivo-only SSDs.
- New tests for missing mouth sizes, reported-first precedence, overrides,
  sampled coefficients (including re-prediction of database-predicted mouth
  sizes, per simulation in `MC_sim_align_parallel()`), `mouth_size_source`,
  seeded SSD fits, the model-domain flag, the intercept standard error,
  in vitro exclusion (`normalize_vivo()`, `ssd_eligible()`, the in vitro
  rainbow trout study, the missing-mouth message) and the corrected bundled
  data.
- `?make_all_pSSDs` and `?run_pssd_reprex` no longer describe `rmore_method =
  "step"` as the method of the published results: the published thresholds
  used `"lognormal"` (Coffin et al. 2026, SI). Defaults are unchanged.
- `data-raw/update_tomex2.R` replaces `data-raw/restore_reported_mouth_size.R`
  and refreshes the bundled data from a ToMEx 2.0 database build.

Manuscript pipeline and repository scripts (outside the package)
----------------------------------------------------------------
The analysis scripts in this repository (`scripts/`), used for Coffin et al.
(2026), received the same fixes:
- `scripts/ToMEx2_EcoTox.Rmd` (deterministic SSDs, quality and sensitivity
  analyses, PSSD++): its inline `align_data()` also defaulted a missing mouth
  size to 5 mm; fixed (44 food dilution records were affected). The threshold
  dataset (`aoc_final`) is restricted to in vivo data (`ssd_eligible()`).
  Its PSSD++ section now defines `erm_registry`, which
  `scripts/utils/functions.R::get_results_df()` has required since it was
  rewritten for the package (commit aaf22cf); without it every PSSD++
  combination failed.
- `scripts/monte carlo/EcoTox_MonteCarlo.Rmd` (MC + SSD): uses the fixed
  `scripts/utils/functions.R::align_data()`, restricts input to in vivo data,
  exports `ssd_eligible()` to the parallel workers, and now sources
  `scripts/utils/plotting_functions.R` (it could not run from a clean session).
- The MC wrappers (`scripts/utils/functions.R::model_wrapper_sobol_parallel()`
  and `MC_sim_align_wrapper()`) read the body-length coefficients as
  `params$sim_beta_log10_body_length` / `params$sim_body_length_intercept`, but
  `matrix_function()` names them `sim.beta.log10.body.length` /
  `sim.body.length.intercept`. The values were `NULL`, so every simulation
  silently used the fixed global coefficients (0.9341, 1.12); the sampled
  coefficients never reached the alignment. Fixed (either spelling is read).
- `model_wrapper_sobol_parallel()` re-aligned data that had already been
  aligned with the default particle-shape parameters. `align_data()` fills a
  missing tested-particle width (length x R.ave) and height (width x H:W
  ratio) only when they are missing, so those stayed at the defaults and the
  sampled R.ave / H:W ratio did not reach tested fragments without a reported
  width (freshwater fragments: up to about 2-fold per simulation). The
  wrapper now aligns the unaligned records each simulation, as
  `PSSDplusplus::MC_sim_align_parallel()` does.
- **Mis-keyed merge of MC summaries (CoV analyses).** Because of the
  re-alignment, the MC results numbered records 1..n of the aligned subset,
  while `aoc_aligned_test` keeps the input row number; `merge(..., by =
  "unique_id")` therefore attached each record's MC mean, SD and coefficient
  of variation to a different record (in the published December 2025 outputs,
  332 of 410 joined records had another record's values; 167 got none). This
  affects the record-level CoV results only: `aoc_MC_sobol.rds`, the CoV
  histograms and `CoV_distribution_summary.csv`, the CoV vs mouth-size tests
  (Spearman, `lm`) and figures (`CoV_mouth_plots`, `CoV_mouth_plotly`), and the
  "MC mean" threshold datasets in `ToMEx2_EcoTox.Rmd`. MC + SSD thresholds and
  PSSD++ inputs use per-simulation values and are not affected. Fixed by the
  change above; the Rmd now stops if the two disagree on any record.
- `scripts/utils/ssd_functions.R`: all SSD filter chains call
  `ssd_eligible()`.
- Checkpointing: the MC Sobol loop saves every 500 simulations to
  `scripts/monte carlo/output/sobol_checkpoints/` and resumes from finished
  chunks (previously a crash lost the whole run of 6-10 h); the PSSD++ loop in
  `ToMEx2_EcoTox.Rmd` reuses a cached combination only if it was computed
  from the same inputs (previously a cache from before an MC re-run would have
  been reused silently); the correction scripts in `docs/correction/R/` save
  each finished step.
- `scripts/utils/alignment_functions.R`: mouth size is predicted only when
  missing (it was recomputed for every record), and `jams_domain` is computed
  with a warning for animals < 30 mm.
- `scripts/utils/ToMEx2.0_Data_Tidying_functions.R` (`ToMEx2.0fxn()`): could
  not run (it used the combined data set before creating it); restructured to
  match the ToMEx 2.0 database script, with the shared
  `scripts/utils/body_length_functions.R` and `scripts/utils/vivo_functions.R`
  (identical to the database repository's copies). An end-to-end test rebuilds
  the bundled data from the ToMEx 2.0 templates and reproduces it exactly.
- `scripts/utils/RDAmaker_functions.R` (ToMEx 1.0 copy): predicted mouth size
  from log10(body length in cm) instead of mm (about 16% too large; not used
  for the bundled data, which come from the database repository) and now
  uses `normalize_vivo()`.
- `scripts/utils/Master_Data_Tidying_Monte_Carlo.R`, `scripts/illustrative_example`:
  5 mm default removed; the former is marked legacy (its source paths no
  longer exist). The illustrative example's text still describes the
  overwritten *Ceriodaphnia dubia* mouth size (104 um; reported 35.7 um) and
  was not re-run.
- `scripts/sync_tomex2_data.R` refuses upstream data that lack the body-length,
  mouth-size and in vivo fields, give a size to unresolved records, or mark
  in vitro records as SSD-eligible (so a sync from a database without
  SCCWRP/ToMEx_AquaticOrganisms#5 cannot silently undo these fixes).
- Script tests in `scripts/tests/`: equivalence of the script and package
  `align_data()` (default and sampled coefficients), body-length functions,
  in vitro exclusion, and the `ToMEx2.0fxn()` rebuild. The deterministic SSD
  inputs of this pipeline were checked against the ToMEx 2.0 Shiny app's SSD
  code (with matching parameters): identical record by record for both ERMs.
- Package vs manuscript pipeline (checked for this release): deterministic
  thresholds for all tiers, environments and ERMs, with their confidence
  limits and species / group / record counts, are identical (88 of 88
  values), whether the package starts from the manuscript's input or from the
  bundled `tomex2`; Monte Carlo aligned doses are identical record by record
  and simulation by simulation; PSSD++ (`step` and `lognormal`) gives
  identical PNEC distributions from identical inputs and seeds. The PSSD++
  plotting grid differs by design (the package spans the observed range).

<!-- Effect on results: to be added after the full re-run (deterministic SSD,
MC + SSD and PSSD++), before the 0.4.0 release is finalized. -->

PSSDplusplus 0.3.5
==================

Bug fixes
---------
- `MC_sim_align_parallel()` now passes `x1D_set` and `x2D_set` through to
  `align_data()` (#5). Previously these bounds were accepted but dropped, so
  every alignment ran at the 1-5000 um defaults regardless of the values
  supplied. Results at the default bounds are unchanged; analyses that set
  other bounds should be re-run.
- New `x1M_set` argument for `MC_sim_align_parallel()` (defaults to
  `x1D_set`) sets the lower bound of the mass/volume/surface-area integrals in
  `align_data()`.
- The internal Sobol wrapper now uses the same size bounds for alignment as
  for unit conversion.
- `align_data()` no longer sets `x2M_ingest` / `x2M_trans` to `NA` for
  records whose `max.size.ingest.um` exactly equals `x2D_set` /
  `upper.tissue.trans.size.um`; the boundary value is now used.

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
