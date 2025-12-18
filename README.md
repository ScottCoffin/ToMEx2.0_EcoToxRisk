# ToMEx 2.0 - Ecotoxicological Risk Assessment

[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D%204.1-blue.svg)](DESCRIPTION)
[![Install PSSDplusplus](https://img.shields.io/badge/install-GitHub-blueviolet.svg)](https://github.com/ScottCoffin/ToMEx2.0_EcoToxRisk/tree/main/package)

Repository for the Journal of Hazardous Materials manuscript ([Pre-print](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5440537)):\
**"A Probabilistic Risk Framework for Microplastics Integrating Uncertainty Across Toxicological and Environmental Variability: Development and Application to Marine and Freshwater Ecosystems."**

Authors: Scott Coffin, Lidwina Bertrand, Kazi Towsif Ahmed, Luan de Souza Leite, Win Cowger, Mariella Sina, Andrew Barrick, Anna Kukkola, Bethanie Carney Almroth, Ezra Miller, Andrew Yeh, Stephanie Kennedy, Magdalena M. Mair\
\
[Related resource: ToMEx Aquatic Organisms](https://github.com/SCCWRP/ToMEx_AquaticOrganisms) (source database and onboarding scripts that underpin the ecotoxicity dataset)

![Graphical Abstract](assets/Graphical_abstract.png)

## Purpose of the assessment

-   Develop PSSD++ thresholds for microplastics, propagating ERM alignment uncertainty, intra-species variability, and parameter uncertainty via Monte Carlo and Sobol sensitivity analysis.
-   Compare aligned toxicity thresholds to environmental particle observations across marine and freshwater compartments.
-   Provide supporting models (tissue translocation GLM) and illustrative materials explaining the alignment framework.

## PSSDplusplus package

To improve the technology readiness level of this method, an R package has been developed to enable rapid adoption with limited barriers. At the moment, we do not intend to host this package on CRAN due to a lack of funding, however we believe that the approach provided here is sufficient for nearly all use cases. If you are interested in helping us adapt this for CRAN, we would be happy to collaborate! 

-   Packaged functions live in `package/` and install as `PSSDplusplus`.
-   Install locally with `devtools::install_local("package", upgrade = "never")`.
-   Load and run the bundled reprex with `library(PSSDplusplus); run_pssd_reprex(sim = 2, n_sim = 2, workers = 1, overwrite_cache = TRUE)`.
-   The walkthrough vignette is available via `vignette("pssdplusplus-walkthrough", package = "PSSDplusplus")`.

Install from GitHub (ignore vignettes by default - as pSSD is computationally heavy):

```r
remotes::install_github("ScottCoffin/ToMEx2.0_EcoToxRisk", subdir = "package", build_vignettes = FALSE)

# or with pak
pak::pak("ScottCoffin/ToMEx2.0_EcoToxRisk", subdir = "package", dependencies = TRUE, build_vignettes =FALSE)
```
The vignette provides a simple walkthrough of how to apply the functions in this package. Due to the high degree of computation required to run even a minimal pSSD++ method, it is recommended to install without building the vignette to avoid unnecessary delays. 

- Load and run the bundled minimal reproducible example:

```r
library(PSSDplusplus)
run_pssd_reprex(sim = 2, n_sim = 2, workers = 1, overwrite_cache = TRUE)
```
- Optional: build/view the vignette after install (can be slow):

```r
remotes::install_github("ScottCoffin/ToMEx2.0_EcoToxRisk", subdir = "package", build_vignettes = TRUE)
vignette("pssdplusplus-walkthrough", package = "PSSDplusplus")
```

### System Requirements

- R 4.1 or newer (matches `Depends` in `DESCRIPTION`).
- Operating systems: standard R platforms (Windows, macOS, Linux)
- Minimum suggested hardware for full Monte Carlo/vignette runs: multi-core CPU and 16 GB RAM; smaller runs (the bundled REPREX) work on modest laptops.
- Positron IDE preferred over RStudio (in general...)

## Repro workflow at a glance

To reproduce the exact analysis in the manuscript:

-   Option A (fastest): Download precomputed Monte Carlo and PSSD++ outputs[from Zenodo](https://doi.org/10.5281/zenodo.16740504) and place them as listed under *Large files*; then knit `scripts/ToMEx2_EcoTox.Rmd`.
-   Option B (full recompute): Run `scripts/monte carlo/EcoTox_MonteCarlo.Rmd` (can exceed 12 hours; high RAM and multiple cores recommended) to generate aligned MC datasets and Sobol outputs, then knit `scripts/ToMEx2_EcoTox.Rmd`.
-   Environmental comparisons: run the scripts in `scripts/characteristics and NMDS/`.
-   Translocation model: knit `scripts/translocation/translocation.Rmd`.
-   Alignment walkthrough: knit `scripts/illustrative_example/ERM Illustrative Example.Rmd`.

## Key scripts and where manuscript components are produced

-   `scripts/ToMEx2_EcoTox.Rmd` — main analysis generating nearly all manuscript figures and tables (threshold tables, PNEC plots, SSD comparisons, etc.). Outputs figures to `output/Manuscript_Figs`, data tables to `output/data`, and threshold/PSDD++ objects to `data/output/`.
-   `scripts/monte carlo/EcoTox_MonteCarlo.Rmd` — Monte Carlo ERM alignment, Sobol sensitivity, aligned dataset creation (sources for MC histograms and threshold distributions). Outputs to `scripts/monte carlo/output/` and selected figures to `output/Manuscript_Figs/`.
-   `scripts/characteristics and NMDS/*.R` — environmental trait extraction and ordinations (ToMEx vs. field particles; Figure 2-type content). Consumes `data/input/Environmental_Data_Collection/00_Extraction_finished/` and writes to `output/Manuscript_Figs/characteristics and NMDS/`.
-   `scripts/translocation/translocation.Rmd` — GLM for tissue translocation probabilities (manuscript translocation section and figures). Outputs knit artifacts and plots to `scripts/translocation/` and `output/Manuscript_Figs/translocation/`.
-   `scripts/illustrative_example/ERM Illustrative Example.Rmd` — didactic alignment example referenced in the text/SI.
-   `scripts/utils/` — shared alignment, SSD, and data-prep functions used by Monte Carlo and PSSD++ workflows.
-   `package/tomex_functions.R` — combined function set (alignment, SSD/PSSD++, tidying) assembled from `scripts/utils/*`.
-   `package/test_tomex_functions.R` — minimal beta test harness using the combined functions (small subset, reduced `nboot`).
-   `assets/` — figures used in the manuscript and supporting visuals.

## Required inputs

-   Toxicity database: `data/input/aoc_z_tomex2.RDS` (originates from the ToMEx repository https://github.com/SCCWRP/ToMEx_AquaticOrganisms and the `aq_mp_tox_shiny` onboarding scripts; see `scripts/monte carlo/readme.txt` for provenance).
-   Legacy/setup data: `data/input/aoc_setup.RDS`, `data/input/aoc_final.RDS` (if present), gape size lists, and related reference files mirrored from ToMEx/`aq_mp_tox_shiny`.
-   Environmental observations: `data/input/Environmental_Data_Collection/00_Extraction_finished/` (used by the characteristics/NMDS scripts).
-   Translocation data: `data/input/translocation_scored_2.xlsx` (and the CSV copy) plus `data/input/aoc_z_tomex2.RDS`.
-   MC reference files: `scripts/monte carlo/ref data/` (e.g., `tomex2_input.rds`, `gape_size.csv`) as noted in `scripts/monte carlo/readme.txt`.

## Dependencies

Install these R packages (R Markdown/knitr required for knitting): - Core analysis (`scripts/ToMEx2_EcoTox.Rmd`): `tidyverse`, `calecopal`, `DT`, `plotly`, `gridExtra`, `grid`, `wesanderson`, `ggtext`, `broom`, `knitr`, `kableExtra`, `viridis`, `ggrepel`, `scales`, `gt`, `ggsci`, `openxlsx`, `ggpubr`, `psych`, `Matrix`, `mc2d`, `trapezoid`, `reshape2`, `devtools`, `sciscales` (install via `devtools::install_github("christyray/sciscales")`), and `ssdtools` version 0.3.7 (checked in-script). - Monte Carlo (`scripts/monte carlo/EcoTox_MonteCarlo.Rmd`): all of the above plus `sensobol`, `truncnorm`, `ggdark`, `gtsummary`, `doParallel`, `doSNOW`, `tictoc`. - Translocation (`scripts/translocation/translocation.Rmd`): `readxl`, `caret`, `DALEX`, `skimr`, `ggeffects`.

## Large files (download from Zenodo or generate locally)

Some outputs are too large for GitHub and are `.gitignore`d. Generate via the scripts or download from `https://doi.org/10.5281/zenodo.16740504` and place with these names: - `scripts/monte carlo/output/`: `sobol_results.rds`, `sobol_results_filtered.rds`, `results_df_sobol.rds`, `all_thresholds_sobol_filtered.rds`, `all_thresholds_sobol_long_filtered.rds`, `mat.rds`, `mat_filtered.rds`, `aoc_MC_sobol.rds`, `summary_stats_base_thresholds_sobol.rds/csv`, `simple_thresholds_table_sobol*.csv`, `simple_RSD.csv`, `sobol_time.txt`. - `data/output/`: `PSDDplusplusresults.rds`, `PNEC_summary_table_wide.rds`, `marine_alpha_thresholds.rds`, tiered SSD files (`tier1_2_*`, `tier3_4_*`). - `output/pssd_cache/` and `output/pssd_debug/`: cached and debug artifacts generated during PSSD++ runs.

## How to reproduce the manuscript results

1)  Install the packages above (confirm `ssdtools` is 0.3.7 and `sciscales` is installed from GitHub).
2)  Place inputs in `data/input/` and MC reference files per `scripts/monte carlo/readme.txt`.
3)  Choose Option A (Zenodo downloads) or Option B (run `scripts/monte carlo/EcoTox_MonteCarlo.Rmd`) to populate MC outputs.
4)  Knit `scripts/ToMEx2_EcoTox.Rmd` to regenerate manuscript figures/tables, PSSD++ results, threshold summaries, and the supplemental workbook (`output/data/supplemental_information.xlsx`).
5)  For environmental comparisons (trait representativeness, NMDS), run `scripts/characteristics and NMDS/`.
6)  For translocation GLM results referenced in the manuscript, knit `scripts/translocation/translocation.Rmd`.
7)  For the pedagogical alignment walkthrough cited in the manuscript/SI, knit `scripts/illustrative_example/ERM Illustrative Example.Rmd`.
8)  To beta-test the packaged functions, run `Rscript package/test_tomex_functions.R` (uses small subset and reduced `nboot`).

## Notes

-   Manuscript text: see the pre-print at https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5440537.
-   The main Rmd (`scripts/ToMEx2_EcoTox.Rmd`) produces the manuscript visuals located in `output/Manuscript_Figs/` (e.g., `Figure5.jpg`, `Figure6_PNEC_compare_arranged_plot.jpg`, `figure1_a_alpha_combined_plot.jpg`, `figure2_bio_response_taxa.jpg`).
-   Environmental data extractions (`data/input/Environmental_Data_Collection/`) underpin comparisons to ToMEx particle traits; keep these synchronized with updates to the ToMEx source data when re-running analyses.

------------------------------------------------------------------------
# ToMEx 2.0 - Ecotoxicological Risk Assessment

Repository for the Journal of Hazardous Materials manuscript ([Pre-print](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5440537)):\
**\"A Probabilistic Risk Framework for Microplastics Integrating Uncertainty Across Toxicological and Environmental Variability: Development and Application to Marine and Freshwater Ecosystems.\"**

Authors: Scott Coffin, Lidwina Bertrand, Kazi Towsif Ahmed, Luan de Souza Leite, Win Cowger, Mariella Sina, Andrew Barrick, Anna Kukkola, Bethanie Carney Almroth, Ezra Miller, Andrew Yeh, Stephanie Kennedy, Magdalena M. Mair\\
\\
![Graphical Abstract](assets/Graphical_abstract.png)

## Purpose of the assessment

-   Develop PSSD++ thresholds for microplastics, propagating ERM alignment uncertainty, intra-species variability, and parameter uncertainty via Monte Carlo and Sobol sensitivity analysis.
-   Compare aligned toxicity thresholds to environmental particle observations across marine and freshwater compartments.
-   Provide supporting models (tissue translocation GLM) and illustrative materials explaining the alignment framework.

## Repro workflow at a glance

-   Option A (fastest): Download precomputed Monte Carlo and PSSD++ outputs[from Zenodo](https://doi.org/10.5281/zenodo.16740504) and place them as listed under *Large files*; then knit `scripts/ToMEx2_EcoTox.Rmd`.
-   Option B (full recompute): Run `scripts/monte carlo/EcoTox_MonteCarlo.Rmd` (can exceed 12 hours; high RAM and multiple cores recommended) to generate aligned MC datasets and Sobol outputs, then knit `scripts/ToMEx2_EcoTox.Rmd`.
-   Environmental comparisons: run the scripts in `scripts/characteristics and NMDS/`.
-   Translocation model: knit `scripts/translocation/translocation.Rmd`.
-   Alignment walkthrough: knit `scripts/illustrative_example/ERM Illustrative Example.Rmd`.

## Key scripts and where manuscript components are produced

-   `scripts/ToMEx2_EcoTox.Rmd` — main analysis generating nearly all manuscript figures and tables (threshold tables, PNEC plots, SSD comparisons, etc.). Outputs figures to `output/Manuscript_Figs`, data tables to `output/data`, and threshold/PSDD++ objects to `data/output/`.
-   `scripts/monte carlo/EcoTox_MonteCarlo.Rmd` — Monte Carlo ERM alignment, Sobol sensitivity, aligned dataset creation (sources for MC histograms and threshold distributions). Outputs to `scripts/monte carlo/output/` and selected figures to `output/Manuscript_Figs/`.
-   `scripts/characteristics and NMDS/*.R` — environmental trait extraction and ordinations (ToMEx vs. field particles; Figure 2-type content). Consumes `data/input/Environmental_Data_Collection/00_Extraction_finished/` and writes to `output/Manuscript_Figs/characteristics and NMDS/`.
-   `scripts/translocation/translocation.Rmd` — GLM for tissue translocation probabilities (manuscript translocation section and figures). Outputs knit artifacts and plots to `scripts/translocation/` and `output/Manuscript_Figs/translocation/`.
-   `scripts/illustrative_example/ERM Illustrative Example.Rmd` — didactic alignment example referenced in the text/SI.
-   `scripts/utils/` — shared alignment, SSD, and data-prep functions used by Monte Carlo and PSSD++ workflows.
-   `package/tomex_functions.R` — combined function set (alignment, SSD/PSSD++, tidying) assembled from `scripts/utils/*`.
-   `package/test_tomex_functions.R` — minimal beta test harness using the combined functions (small subset, reduced `nboot`).
-   `assets/` — figures used in the manuscript and supporting visuals.

## Required inputs

-   Toxicity database: `data/input/aoc_z_tomex2.RDS` (originates from the ToMEx repository https://github.com/SCCWRP/ToMEx_AquaticOrganisms and the `aq_mp_tox_shiny` onboarding scripts; see `scripts/monte carlo/readme.txt` for provenance).
-   Legacy/setup data: `data/input/aoc_setup.RDS`, `data/input/aoc_final.RDS` (if present), gape size lists, and related reference files mirrored from ToMEx/`aq_mp_tox_shiny`.
-   Environmental observations: `data/input/Environmental_Data_Collection/00_Extraction_finished/` (used by the characteristics/NMDS scripts).
-   Translocation data: `data/input/translocation_scored_2.xlsx` (and the CSV copy) plus `data/input/aoc_z_tomex2.RDS`.
-   MC reference files: `scripts/monte carlo/ref data/` (e.g., `tomex2_input.rds`, `gape_size.csv`) as noted in `scripts/monte carlo/readme.txt`.

## Dependencies

Install these R packages (R Markdown/knitr required for knitting): - Core analysis (`scripts/ToMEx2_EcoTox.Rmd`): `tidyverse`, `calecopal`, `DT`, `plotly`, `gridExtra`, `grid`, `wesanderson`, `ggtext`, `broom`, `knitr`, `kableExtra`, `viridis`, `ggrepel`, `scales`, `gt`, `ggsci`, `openxlsx`, `ggpubr`, `psych`, `Matrix`, `mc2d`, `trapezoid`, `reshape2`, `devtools`, `sciscales` (install via `devtools::install_github(\"christyray/sciscales\")`), and `ssdtools` version 0.3.7 (checked in-script). - Monte Carlo (`scripts/monte carlo/EcoTox_MonteCarlo.Rmd`): all of the above plus `sensobol`, `truncnorm`, `ggdark`, `gtsummary`, `doParallel`, `doSNOW`, `tictoc`. - Translocation (`scripts/translocation/translocation.Rmd`): `readxl`, `caret`, `DALEX`, `skimr`, `ggeffects`.

## Large files (download from Zenodo or generate locally)

Some outputs are too large for GitHub and are `.gitignore`d. Generate via the scripts or download from `https://doi.org/10.5281/zenodo.16740504` and place with these names: - `scripts/monte carlo/output/`: `sobol_results.rds`, `sobol_results_filtered.rds`, `results_df_sobol.rds`, `all_thresholds_sobol_filtered.rds`, `all_thresholds_sobol_long_filtered.rds`, `mat.rds`, `mat_filtered.rds`, `aoc_MC_sobol.rds`, `summary_stats_base_thresholds_sobol.rds/csv`, `simple_thresholds_table_sobol*.csv`, `simple_RSD.csv`, `sobol_time.txt`. - `data/output/`: `PSDDplusplusresults.rds`, `PNEC_summary_table_wide.rds`, `marine_alpha_thresholds.rds`, tiered SSD files (`tier1_2_*`, `tier3_4_*`). - `output/pssd_cache/` and `output/pssd_debug/`: cached and debug artifacts generated during PSSD++ runs.

## PSSD++ minimal REPREX (package/PSSD_reprex.R)

A compact, end-to-end reproducible example that runs the PSSD++ workflow used in the manuscript. The script lives at `package/PSSD_reprex.R` and demonstrates the full Monte Carlo alignment -> pSSD pipeline using reduced simulation sizes for fast iteration and debugging.

Key points

- The script assembles required helper functions from `scripts/utils/` and `scripts/PSSD/` and enforces `ssdtools` v0.3.7 for consistent SSD behaviour.
- Defaults are intentionally small for the REPREX (`n_sobol = 10`, Monte Carlo `n_sim = 10`, PSSD `sim = 10`). Increase these (e.g. `sim >= 300`) for production-quality results.
- Results are cached and figures are saved so repeated runs avoid recomputation (see cache and output paths below).

Inputs required

- `data/input/aoc_z_tomex2.RDS` — preprocessed toxicity dataset used by the example.
- Helper scripts in `scripts/utils/` and `scripts/PSSD/` (these are sourced by the REPREX script).

Where outputs go

- Cache: `package/test_output/pssd_cache/` (per combination .rds cache files)
- Figures: `package/test_output/figures/` (publication-ready PSSD and PNEC plots)

How to run

Run non-interactively from a shell:

```bash
Rscript package/PSSD_reprex.R
```

Or source interactively from an R session:

```r
source("package/PSSD_reprex.R")
```

Notes & tips

- The script will check and attempt to install the required `ssdtools` version; ensure you have permissions to install packages if needed.
- To increase fidelity change the `sim`, `n_sim` (Monte Carlo sample size) and `workers` variables inside the script before running.
- Set `overwrite_cache = TRUE` in the `make_all_pSSDs()` call to force recomputation of cached results.

Expected output

- Cached pSSD `.rds` files and saved figures as described above.
- The script prints a summarized PNEC table (`PNEC_summary`) at the end and leaves full pSSD objects in the cache for downstream analysis or plotting.

## How to reproduce the manuscript results

1)  Install the packages above (confirm `ssdtools` is 0.3.7 and `sciscales` is installed from GitHub).
2)  Place inputs in `data/input/` and MC reference files per `scripts/monte carlo/readme.txt`.
3)  Choose Option A (Zenodo downloads) or Option B (run `scripts/monte carlo/EcoTox_MonteCarlo.Rmd`) to populate MC outputs.
4)  Knit `scripts/ToMEx2_EcoTox.Rmd` to regenerate manuscript figures/tables, PSSD++ results, threshold summaries, and the supplemental workbook (`output/data/supplemental_information.xlsx`).
5)  For environmental comparisons (trait representativeness, NMDS), run `scripts/characteristics and NMDS/`.
6)  For translocation GLM results referenced in the manuscript, knit `scripts/translocation/translocation.Rmd`.
7)  For the pedagogical alignment walkthrough cited in the manuscript/SI, knit `scripts/illustrative_example/ERM Illustrative Example.Rmd`.
8)  To beta-test the packaged functions, run `Rscript package/test_tomex_functions.R` (uses small subset and reduced `nboot`).

## Notes

-   Manuscript text: see the pre-print at https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5440537.
-   The main Rmd (`scripts/ToMEx2_EcoTox.Rmd`) produces the manuscript visuals located in `output/Manuscript_Figs/` (e.g., `Figure5.jpg`, `Figure6_PNEC_compare_arranged_plot.jpg`, `figure1_a_alpha_combined_plot.jpg`, `figure2_bio_response_taxa.jpg`).
-   Environmental data extractions (`data/input/Environmental_Data_Collection/`) underpin comparisons to ToMEx particle traits; keep these synchronized with updates to the ToMEx source data when re-running analyses.

------------------------------------------------------------------------
