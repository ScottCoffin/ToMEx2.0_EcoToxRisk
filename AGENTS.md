# AGENTS.md — working in ToMEx2.0_EcoToxRisk

Guidance for people and coding agents changing this repository. Read this
before editing any alignment, SSD or PSSD++ code. Keep it current when you
learn something the hard way.

## What is here

Two implementations of the same science live side by side:

| Code base | Used for | Entry points |
|---|---|---|
| `scripts/` (sourced R files + Rmd) | The published analysis, Coffin et al. 2026 (J. Hazard. Mater.) | `scripts/ToMEx2_EcoTox.Rmd`, `scripts/monte carlo/EcoTox_MonteCarlo.Rmd` |
| `package/` (R package **PSSDplusplus**) | Reuse by others (e.g. `../W2W_river_plastics`) | `align_data()`, `MC_sim_align_parallel()`, `make_all_pSSDs()`, `make_tiered_SSDs()`, `hc5_from_endpoints()` |

**They must give identical numbers.** A fix in one is a fix in both, followed
by the equivalence tests (below). The upstream database code
(`SCCWRP/ToMEx_AquaticOrganisms`, the Shiny app) shares
`scripts/utils/body_length_functions.R`, `vivo_functions.R` and the
deterministic alignment; keep those copies identical too.

## Data flow (sitemap)

```
ToMEx 2.0 database (SCCWRP/ToMEx_AquaticOrganisms, tidying + body-length/vivo corrections)
  └─ scripts/sync_tomex2_data.R  (refuses data missing the v0.4.0 fields)
       └─ data/input/aoc_z_tomex2.RDS ──┬─ package/data-raw/update_tomex2.R → package/data/tomex2.rda (+ tomex2_mini, inst/extdata)
                                        ├─ scripts/monte carlo/EcoTox_MonteCarlo.Rmd   (MC alignment + per-sim SSDs; ~6.5 h on 12 cores)
                                        │     └─ scripts/monte carlo/output/{mat,sobol_results,results_df_sobol,aoc_MC_sobol}.rds
                                        └─ scripts/ToMEx2_EcoTox.Rmd  (deterministic SSDs, quality/sensitivity, CoV, PSSD++)
                                              reads results_df_sobol.rds (PSSD++ input) and aoc_MC_sobol.rds (CoV analyses)
                                              PSSD++ cache: output/pssd_cache/<rmore_method>_<sim>sims/<combo>.rds
                                              figures/tables: output/Manuscript_Figs, output/data, data/output
```

Folders not on the main path: `scripts/translocation/` (tissue translocation
GLM → `beta_0`, `beta_1`), `scripts/characteristics and NMDS/`,
`scripts/illustrative_example/` (text not re-run after 0.4.0), and
`scripts/utils/Master_Data_Tidying_Monte_Carlo.R` (legacy; its paths no longer
exist).

## Mirrored functions (edit both)

| Package (`package/R/`) | Script (`scripts/utils/`) |
|---|---|
| `functions.R::align_data()` | `functions.R::align_data()` (also an inline copy in `ToMEx2_EcoTox.Rmd`) |
| `functions.R::matrix_function()` | `functions.R::matrix_function()` (**defined twice**; the second, ~line 1270, wins) |
| `functions.R::MC_sim_align_wrapper()` / `MC_sim_align_parallel()` | `functions.R` same names; MC Rmd uses `model_wrapper_sobol_parallel()` |
| `ssd_functions.R::SSD_function_t*()`, `process_environment_data()` | `ssd_functions.R` same names |
| `functions.R::generate_plots_and_summary()`, `prep_data()`, `prepare_plot_data()` | `functions.R` same names; `scripts/PSSD/{rmore,do.pssd}.r` |
| `vivo.R::ssd_eligible()`, `normalize_vivo()` | `vivo_functions.R` |

Known, intended difference: the package's PSSD++ plotting grid
(`$quantiles$iv`) spans the observed range instead of 10^-5–10^10 (plots
only). Both seed each PSSD++ combination with `2026 + sum(utf8ToInt(combo_id))`
(`make_all_pSSDs(combo_seed = TRUE)`, the default), so they give identical
draws whatever the order of combinations or the caller's seed.

## Verify a change

```sh
# package tests (renv is active inside package/; run from there)
cd package && Rscript -e 'pkgload::load_all(); testthat::test_dir("tests/testthat")'
# script tests (repo root)
Rscript -e "testthat::test_file('scripts/tests/test_align_data_equivalence.R')"   # script vs package align_data
Rscript -e "testthat::test_file('scripts/tests/test_vivo_exclusion.R')"
Rscript -e "testthat::test_file('scripts/tests/test_body_length_functions.R')"
TOMEX2_RAW_RDS=<raw build> Rscript -e "testthat::test_file('scripts/tests/test_tomex2_tidying.R')"
```

After changing alignment/SSD/PSSD++ code, also compare package vs Rmd on real
data. The 0.4.0 check (deterministic thresholds with CIs; MC aligned doses per
record and simulation; PSSD++ PNECs from identical inputs and seed) found
identical results; repeat it if you touch those paths.

## Traps that have cost real time

- **Hidden globals in `scripts/`.** Script functions read objects the Rmd
  defines: `nboot`, `x1M_set`, `x1D_set`, `x2D_set`, `n_sim`,
  `sim.upper.tissue.trans.size.um` (set inside the MC `foreach` body),
  `aoc_risk_paper`, `aoc_aligned`, `beta_log10_body_length`,
  `global_color_palette`, `erm_registry` (read by `get_results_df()`). A
  missing global errors; a stale one silently changes results. `foreach`
  only sees them if listed in `.export`.
- **Render the Rmd in the global environment** (`rmarkdown::render(...,
  envir = globalenv())`); the leave-one-out code relies on the sourced
  `utils/ssd_functions.R` replacing the Rmd's own `SSD_function_t*`.
- **ssdtools must be 0.3.7** to reproduce the manuscript (2.x gives different
  deterministic HCx). With 0.3.7, ssdtools must also be *attached*: fitdistrplus
  finds distributions on the search path, and an unattached ssdtools silently
  drops them from the model average. The package now attaches it; any new
  script code must `library(ssdtools)`. PSSD++ does not use ssdtools.
- **Parallel workers are fresh R sessions** and do not inherit library paths
  set with `.libPaths()` inside the main session: they load whatever is in
  the default library (e.g. ssdtools 2.5.0), silently changing every SSD.
  The MC Rmd and the package now pass `.libPaths()` to their workers, and the
  MC Rmd stops unless every worker has ssdtools 0.3.7. Setting `R_LIBS`
  before launching R also works.
- **Seeds.** SSD fits use `set.seed(99)` (Tier 1 = HC5 lower confidence limit
  from only `nboot = 10` bootstraps, so it is seed-sensitive). Each PSSD++
  combination is seeded from its name (Rmd loop and `make_all_pSSDs()`).
- **`unique_id`** is the row number of the data *passed to* `align_data()`,
  assigned before its internal filtering. Re-aligning already-aligned data
  renumbers rows (and freezes widths/heights filled from the default R.ave and
  H:W ratio). Always align unaligned records, and join MC output to records by
  `unique_id` only when both sides came from the same input.
- **Mouth size precedence:** reported (record) > reported (species, life
  stage) > `max.size.ingest.override.um` > predicted from body length (Jams et
  al. 2020) > unresolved (excluded from Food Dilution; never 5 mm). In
  `tomex2`, `max.size.ingest.um` also stores the database's own predictions;
  use `max.size.ingest.basis` to tell reported from predicted.
- **`align_data()` drops rows** (12,797 → ~9,500 for `tomex2`); match on
  `rowid`, not position.
- **renv** is active only inside `package/`; the user library outside it has
  different versions (e.g. ssdtools 2.5.0). The Rmd's `check_pckg.R` may try
  to install pinned versions.
- **Runtimes:** MC Sobol ~6.5 h on 12 cores. PSSD++ with `rmore_method =
  "lognormal"` (the published method, see the SI) takes ~40 min per
  combination at `sim = 300`; the trapezoidal `"step"` method took >8.5 h for
  one combination. A PSSD++ run at production size (9,424 MC
  simulations x `sim = 300`) commits ~60 GB: run it **alone** on a 64 GB
  machine (not alongside the MC, a render or another PSSD++ run), or the
  system runs out of commit memory and other jobs crash. Run combinations as separate
  processes (memory, not cores, is the limit: ~64 GB for 300 sims).
- **Render the PSSD++ loop in separate processes.** In one `rmarkdown::render`
  session the loop's memory grew to ~180 GB committed by the fourth
  combination (2026-10-21) and the machine started paging heavily. Run the
  Rmd up to the loop, save the workspace, run each combination in its own R
  process with `pssd_only <- "<combo_id>"` (the loop skips the others), then
  render: the loop loads the cached combinations (same inputs, same
  fingerprint). Each combination has its own seed (`2026 + sum(utf8ToInt(combo_id))`),
  so results are identical either way.
- **Checkpointing (resume after a crash: re-run the same command).** MC Sobol
  loop: `scripts/monte carlo/output/sobol_checkpoints/` (every 500 sims);
  manuscript PSSD++: `output/pssd_cache/<method>_<sim>sims/` (per
  combination); package: `MC_sim_align_parallel(checkpoint_dir =)` and
  `make_all_pSSDs(base_cache_dir =)`; correction attribution:
  `docs/correction/results/checkpoints/`. All are keyed by a fingerprint of
  their inputs, so stale results are recomputed, never reused. Delete a
  checkpoint folder to force a clean run. Inside `foreach`, keep the call at
  top level (a for loop, not `lapply(function ...)`), or the automatic export
  of helper functions breaks.
- **Never edit a script while Rscript or bash is running it**: both read the
  file as they go, so an edit corrupts the running job. Stage changes in a
  copy and swap them in when the run has finished.
- **Large outputs live on OneDrive** (`output/` ~1 GB). Keep caches and
  scratch runs out of git (`output/pssd_cache*`, `scripts/.ToMEx2_EcoTox_rerun.Rmd`).
- Writing R files that contain backslashes through a shell heredoc can strip
  them; use an editor/file-write tool.

## To-do: production runs and release (keep this list current)

Order matters; each step needs the previous one.

- [x] **Re-run the MC**: done 2026-10-21 12:03 (worker check passed; an
      earlier 06:27 run was invalid, its workers loaded ssdtools 2.5.0).
- [x] `nboot` decision: stays 10 (analysis, package, Shiny app); the effect
      of 1000 is tested on a subsample by `docs/correction/R/nboot_impact.R`
      and reported in the correction (Table C4).
- [x] **Re-render `scripts/ToMEx2_EcoTox.Rmd`**: done 2026-10-21 15:23
      (PSSD++ combinations in separate processes, then the render from the
      cache; no chunk errors). A `"step"` PSSD++ run (~90 h) is optional.
- [x] **Attribution**: done 2026-10-22 with 1,000 MC x 300 PSSD++ sims per
      step (decided with the author), steps run one at a time. Full-size
      (9,424) S0-S2 agree within 4% (`docs/correction/results/fullsize_check/`);
      S6 reproduces the production render within 9%. Correction rendered.
- [x] `docs/correction/R/cov_misjoin_impact.R` on the published MC output.
- [x] CoV statistics of the corrected MC added to the correction.
- [ ] **Archive the published outputs** that are not in git
      (`scripts/monte carlo/output/pre_issue6/`, `output/pssd_cache_pre_issue6/`)
      to a permanent repository (e.g. Zenodo) and cite them in the correction.
      Needs the author's account.
- [x] Draft "Effect on the conclusions" (correction) and "Effect on results"
      (`package/NEWS.md`) from the production render. **Authors to review.**
- [x] **Release 0.4.0**: code (58114b0) plus regenerated outputs and the
      final correction committed on branch `v0.4.0-corrections` and pushed
      2026-10-22. Next: review/merge into `main` (closes #6); authors to review
      the drafted conclusions; archive the published outputs (item above).

## Open technical debt

- The package's internal `model_wrapper_sobol_parallel()` (MC + per-sim SSDs)
  is unused and does not run: it reads `aoc_aligned`, `x2D_set` and
  `sim.upper.tissue.trans.size.um` from the global environment and then fails
  in `convert_units_fxn()`. Either make it a working exported function with
  explicit arguments or delete it.
- Duplicate `matrix_function()` definitions in `scripts/utils/functions.R`.
- `scripts/` functions should take explicit arguments instead of globals;
  converting the script pipeline to call the package would remove the mirror
  entirely.
- 0.4.0 release: see the to-do list above.
