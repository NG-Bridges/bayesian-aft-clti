#' Analysis settings shared by every script.
#' Values that enter `spec_string()` key the saved fits; changing one means refitting everything.
analysis_config <- function() {
  cfg <- list(
    seed = 2026,
    horizons_months = c(3, 6, 12, 24, 36),
    n_folds = 5,
    n_imputations = 50, # m >= 100 x FMI, with FMI 0.50 for GNRI (von Hippel 2020)
    rsf_n_imputations = 50,
    cv_n_imputations = 50,
    # Cross-validating the frailty and mixture models would add ~1,000 fits (weeks);
    # they are compared with the single-component model on WAIC/DIC3 instead.
    cv_models = "AFT_single",
    models = c("AFT_single", "AFT_mixture", "AFT_frailty"),
    stan_chains = 4, stan_cores = 4,
    stan_warmup = 1500, stan_sample = 1500,
    stan_adapt_delta = 0.99, stan_treedepth = 12,
    intercept_scale = 5,
    shape_scale = 0.5,
    design_scale = 2,
    gate_scale = 2,
    slab_scale = 1, slab_df = 4,
    p0_relevant = 5,
    frailty_scale = 1,
    # WAIC within ~0.7 of a converged 100-node value at the amputation posterior (sigma_u 2-3.6),
    # though single draws can be off by several deviance units; checked in Phase 1
    gh_nodes = 25,
    K_mix = 2,
    relabel_draws = 2000, # Stephens holds relabel_draws x N x K at once (~54 MB)
    anchor_var = "eGFR",
    include_cohort = TRUE,
    include_egfr_sq = TRUE,
    # 1,500 draws: |dAUC| ~ 1e-3 against a 6,000-draw reference, at a few % of Phase 3 runtime
    pred_draws = 1500,
    grid_n = 150,
    store_draws = 1000,
    rsf_ntree = 1000, rsf_nodesize = 15, rsf_nsplit = 10,
    rsf_cores = max(1, parallel::detectCores() - 1),
    rhat_max = 1.01, ess_min = 400, div_frac_max = 0.005,
    calib_models = c("AFT_single", "RSF"),
    calib_n_groups = 10,
    calib_spline_df = 3,
    calib_n_boot = 200,
    calib_main_months = c(12, 36),
    calib_verify_n = 3,
    ari_prob_cut = 0.80,
    # A partition read off chains that have not mixed is not evidence; the amputation
    # mixture reaches rhat_lp < 1.01 in about half its fits, so ARI is also reported on that subset
    ari_rhat_lp_max = 1.01,
    leak_boot_B = 200,
    leak_boot_months = 24
  )
  if (example_mode()) cfg <- utils::modifyList(cfg, example_settings())
  cfg$horizons <- cfg$horizons_months / 12
  cfg
}

#' TRUE when CLTI_EXAMPLE=1: run on a simulated export instead of the trial data.
example_mode <- function() identical(Sys.getenv("CLTI_EXAMPLE"), "1")

#' Smaller settings for the example run; its fits are keyed and stored apart from the real ones.
example_settings <- function() {
  list(
    n_imputations = 3, rsf_n_imputations = 3, cv_n_imputations = 3,
    stan_chains = 2, stan_cores = 2, stan_warmup = 300, stan_sample = 300,
    store_draws = 300, pred_draws = 300, relabel_draws = 300,
    rsf_ntree = 200, calib_n_boot = 20, leak_boot_B = 20
  )
}

#' Where the pipeline reads and writes. The example run never touches the real data, fits or outputs.
project_paths <- function() {
  if (example_mode()) {
    list(
      derived = here::here("data", "derived", "example"),
      results = here::here("results", "example"), output = here::here("output", "example")
    )
  } else {
    list(
      raw = here::here("data", "raw"), derived = here::here("data", "derived"),
      results = here::here("results"), output = here::here("output")
    )
  }
}

#' Path inside the results directory.
results_path <- function(...) file.path(project_paths()$results, ...)

#' Competing causes: event column and status code.
cause_codes <- function() {
  list(
    amputation = list(event = "ev_amp", code = 1L),
    death = list(event = "ev_death", code = 2L)
  )
}
