#' Run the PSSD++ reproducible example
#'
#' Executes the full Monte Carlo alignment and PSSD++ workflow bundled with the
#' repository, using the included ToMEx 2.0 ecotoxicity dataset.
#'
#' @param data_path Path to the input dataset. Defaults to the bundled
#'   `aoc_z_tomex2.RDS`.
#' @param cache_dir Directory for cached PSSD objects.
#' @param figure_dir Directory for saved figure outputs.
#' @param sim Number of PSSD simulations per combination.
#' @param n_sim Number of Monte Carlo alignment simulations.
#' @param workers Number of workers for both alignment and PSSD phases.
#' @param overwrite_cache Logical; force recalculation even if cache exists.
#' @param parallel Logical; run PSSD combinations in parallel.
#' @param seed Optional integer seed for reproducibility.
#' @param quantile_type Quantile type used when summarizing outputs.
#' @param cv_uf Coefficient of variation for uncertainty factors.
#' @param rmore_method Distribution method for sampling species-level NOEC
#'   distributions in the PSSD++ step. `"step"` (default) uses a trapezoidal
#'   piecewise-uniform approximation following Wigger et al. (2020), which is
#'   consistent with published results. `"lognormal"` uses a log-normal
#'   short-cut that is faster but assumes symmetric log-scale variability; use
#'   only when speed matters more than strict comparability to the manuscript.
#'
#' @return A list containing the aligned data, PSSD outputs, and summary table.
#' @export
#'
#' @examples
#' \dontrun{
#'   # Fastest possible run using the bundled mini dataset (seconds)
#'   run_pssd_reprex(sim = 1, n_sim = 1, workers = 1, overwrite_cache = TRUE)
#'
#'   # Reproducible run with a fixed seed
#'   result <- run_pssd_reprex(sim = 2, n_sim = 2, workers = 1, seed = 42)
#'   result$pnec_summary
#'
#'   # Use lognormal method (faster, less precise)
#'   run_pssd_reprex(sim = 2, n_sim = 2, workers = 1, rmore_method = "lognormal")
#' }
run_pssd_reprex <- function(
  data_path = NULL,
  cache_dir = file.path(tempdir(), "pssd_cache"),
  figure_dir = file.path(tempdir(), "pssd_figures"),
  sim = 10,
  n_sim = 10,
  workers = max(1, parallel::detectCores() - 1),
  overwrite_cache = FALSE,
  parallel = TRUE,
  seed = NULL,
  quantile_type = 8,
  cv_uf = 0.5,
  rmore_method = "step"
) {
  # Load bundled dataset if no custom path is provided
  if (is.null(data_path)) {
    data_path <- system.file("extdata", "aoc_z_tomex2.RDS", package = "PSSDplusplus")
  }

  if (!nzchar(data_path) || !file.exists(data_path)) {
    # Fall back to the attached data object if the file is missing
    if (exists("tomex2", envir = asNamespace("PSSDplusplus"))) {
      aoc_risk_paper <- get("tomex2", envir = asNamespace("PSSDplusplus"))
    } else {
      stop("Input data file not found. Specify a valid `data_path`.", call. = FALSE)
    }
  }

  check_ssdtools_version()
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Prepare toxicity dataset
  if (!exists("aoc_risk_paper")) {
    aoc_risk_paper <- readRDS(data_path)
  }

  aoc_risk_paper <- aoc_risk_paper %>%
    dplyr::rename(environment = env_f) %>%
    dplyr::mutate(AF.total = af.time * af.noec) %>%
    tidyr::drop_na(effect.metric) %>%
    dplyr::filter(
      !environment %in% c("Terrestrial", "Not Reported"),
      Group != "Bacterium",
      Group != "Plant",
      effect.metric != "HONEC",
      tier_zero_tech_f == "Red Criteria Passed",
      tier_zero_risk_f == "Red Criteria Passed",
      risk.13 != 0
    )

  # Generate parameter matrix
  mat <- matrix_function(
    n_sobol = n_sim,
    params = param_default_values,
    upper.tissue.truncation.limit = 500,
    x1M_set = 1,
    x2D_set = 5000
  )

  parameter_histograms <- parameter_histograms_function(mat)

  # Alignment simulations
  MC_sim_df <- MC_sim_align_parallel(
    tox_data = aoc_risk_paper,
    param_matrix = mat,
    x1D_set = 1,
    x2D_set = 5000,
    n_sim = n_sim,
    num_cores = workers
  )

  # Build ERM-specific datasets
  results_df_food <- MC_sim_df |>
    dplyr::filter(
      ingestible != "not ingestible",
      particles_L_food_dilution > 0,
      Group != "Algae"
    ) |>
    dplyr::mutate(dose_new_particles_L = particles_L_food_dilution) |>
    tidyr::drop_na(particles_L_food_dilution)
  results_df_food_t3_t4 <- results_df_food |>
    dplyr::filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))

  results_df_tissue <- MC_sim_df |>
    dplyr::filter(translocatable != "not translocatable", particles_L_ox_stress > 0) |>
    dplyr::mutate(dose_new_particles_L = particles_L_ox_stress) |>
    tidyr::drop_na(particles_L_ox_stress)
  results_df_tissue_t3_t4 <- results_df_tissue |>
    dplyr::filter(risk.13 != 1, bio_f %in% c("Organism", "Population"))

  erm_registry <- list(
    "Food Dilution" = list(
      base = results_df_food,
      t3_t4 = results_df_food_t3_t4
    ),
    "Tissue Translocation" = list(
      base = results_df_tissue,
      t3_t4 = results_df_tissue_t3_t4
    )
  )

  pSSDs <- make_all_pSSDs(
    MC_sim_df = MC_sim_df,
    tiers = c(3),
    environments = c("Freshwater", "Marine"),
    erms = c("Food Dilution", "Tissue Translocation"),
    erm_registry = erm_registry,
    sim = sim,
    cv_uf = cv_uf,
    rmore_method = rmore_method,
    quantile_type = quantile_type,
    debug_option = FALSE,
    parallel = parallel,
    workers = workers,
    base_cache_dir = cache_dir,
    base_output_path = figure_dir,
    overwrite_cache = overwrite_cache
  )

  PNEC_summary <- summarize_PNECs(pSSDs)

  list(
    data = aoc_risk_paper,
    param_matrix = mat,
    parameter_histograms = parameter_histograms,
    MC_sim_df = MC_sim_df,
    erm_registry = erm_registry,
    pssd_results = pSSDs,
    pnec_summary = PNEC_summary
  )
}
