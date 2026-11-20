#' Read the prepared data written by `01_prepare_data.R`.
#' @return List: `P` (predictor sets), `base` (one row per patient), `imps` (completed datasets).
load_prepped <- function(path = file.path(project_paths()$derived, "prepped.rds")) {
  if (!file.exists(path)) stop("Missing ", path, "; run scripts/01_prepare_data.R first.")
  s1 <- readRDS(path)
  need <- c("cr_time", "cr_status", "ev_amp", "ev_death")
  if (!all(need %in% names(s1$base_ds))) {
    stop("prepped data lacks: ", paste(setdiff(need, names(s1$base_ds)), collapse = ", "))
  }
  list(P = s1$P, base = s1$base_ds, imps = s1$imputed)
}

#' Split predictors into the unpenalized design block, the horseshoe block and the mixture gate.
predictor_sets <- function(core, cfg) {
  design_vars <- c("ARM", "FACTOR_B", if (cfg$include_cohort) "COHORT")
  pen_vars <- setdiff(unique(c(core, "GNRI")), design_vars)
  model_pen_vars <- if (cfg$include_egfr_sq) c(pen_vars, "eGFR_sq") else pen_vars
  gate_vars <- intersect(c("eGFR", "AGE_T", "AGE", "AMB_STAT3CAT", "DIAB_YN", "GNRI"), pen_vars)
  if ("AGE_T" %in% gate_vars) gate_vars <- setdiff(gate_vars, "AGE")
  list(
    design_vars = design_vars, pen_vars = pen_vars,
    model_pen_vars = model_pen_vars, gate_vars = gate_vars
  )
}

#' Check that every predictor is present, numeric and complete before hours of fitting.
check_predictors <- function(vs, base, imps) {
  all_vars <- c(vs$pen_vars, vs$design_vars)
  if (length(miss <- setdiff(all_vars, names(base)))) {
    stop("Predictors absent from the data: ", paste(miss, collapse = ", "))
  }
  numeric_ok <- vapply(base[, all_vars, drop = FALSE], function(z) is.numeric(z) || is.logical(z), logical(1))
  if (!all(numeric_ok)) stop("Non-numeric predictors: ", paste(all_vars[!numeric_ok], collapse = ", "))
  # Checked on a completed dataset: base still carries raw albumin, so GNRI is partly missing there
  na_frac <- vapply(imps[[1]][, all_vars, drop = FALSE], function(z) mean(is.na(z)), numeric(1))
  if (any(na_frac > 0)) stop("Missing values in completed dataset 1: ", paste(all_vars[na_frac > 0], collapse = ", "))
  if (length(vs$gate_vars) < 2) stop("Fewer than 2 gate covariates available.")
  invisible(TRUE)
}

#' Recode binary design variables to 0/1, the higher raw level becoming 1.
recode_design <- function(dfs, vars, ref) {
  for (v in vars) {
    lv <- sort(unique(ref[[v]]))
    if (length(lv) != 2) stop(v, " has ", length(lv), " levels; design variables must be binary.")
    for (i in seq_along(dfs)) dfs[[i]][[v]] <- as.integer(dfs[[i]][[v]] == lv[2])
    ref[[v]] <- as.integer(ref[[v]] == lv[2])
  }
  list(dfs = dfs, ref = ref)
}

std_moments <- function(df, vars) {
  list(mu = sapply(df[, vars, drop = FALSE], mean), sd = sapply(df[, vars, drop = FALSE], sd))
}

std_apply <- function(df, vars, mom) {
  df[, vars] <- sweep(sweep(df[, vars, drop = FALSE], 2, mom$mu, "-"), 2, mom$sd, "/")
  df
}

#' Fixed standardization moments from the stacked completed datasets (inference fits).
fixed_moments <- function(imps, vars) {
  stacked <- do.call(rbind, lapply(imps, function(d) d[, vars, drop = FALSE]))
  mom <- std_moments(stacked, vars)
  if (any(mom$sd <= 0)) stop("Zero-variance predictors: ", paste(vars[mom$sd <= 0], collapse = ", "))
  mom
}

#' Standardize, then add the square of standardized eGFR.
prepare_design <- function(df, vs, mom, cfg) {
  df <- std_apply(df, vs$pen_vars, mom)
  if (cfg$include_egfr_sq) df$eGFR_sq <- df$eGFR^2
  df
}

design_mat <- function(df, vars) {
  m <- as.matrix(df[, vars, drop = FALSE])
  storage.mode(m) <- "double"
  m
}

#' One fold partition, stratified on competing-risk status, shared by every model and dataset.
make_folds <- function(strata, k, seed) {
  set.seed(seed)
  fold <- integer(length(strata))
  for (lv in unique(strata)) {
    idx <- sample(which(strata == lv))
    fold[idx] <- rep_len(seq_len(k), length(idx))
  }
  fold
}

#' Global horseshoe scale per cause (Piironen and Vehtari 2017).
#' sigma comes from an intercept-only Weibull fit; Gumbel scale 1/k has SD pi / (k sqrt(6)).
horseshoe_tau0 <- function(base, n_pen, cfg) {
  vapply(cause_codes(), function(cs) {
    f0 <- survival::survreg(survival::Surv(base$cr_time, base[[cs$event]]) ~ 1, dist = "weibull")
    sig <- f0$scale * pi / sqrt(6)
    p0 <- min(cfg$p0_relevant, n_pen - 1)
    p0 / (n_pen - p0) * sig / sqrt(nrow(base))
  }, numeric(1))
}

#' Gauss-Hermite nodes and weights (Golub-Welsch), physicists' convention.
gauss_hermite <- function(n) {
  i <- seq_len(n - 1)
  b <- sqrt(i / 2)
  jac <- matrix(0, n, n)
  jac[cbind(i, i + 1)] <- b
  jac[cbind(i + 1, i)] <- b
  e <- eigen(jac, symmetric = TRUE)
  ord <- order(e$values)
  gh <- list(x = e$values[ord], w = sqrt(pi) * (e$vectors[1, ord])^2)
  gh$logw <- log(gh$w / sqrt(pi))
  stopifnot(abs(sum(exp(gh$logw)) - 1) < 1e-8)
  gh
}

#' Everything a fit or prediction needs, built once: predictor sets, moments, tau0, quadrature.
model_setup <- function(dat, cfg) {
  vs <- predictor_sets(dat$P$core, cfg)
  check_predictors(vs, dat$base, dat$imps)
  rc <- recode_design(dat$imps, vs$design_vars, dat$base)
  c(vs, list(
    cfg = cfg,
    base = rc$ref,
    imps = rc$dfs,
    fix_mom = fixed_moments(rc$dfs[seq_len(cfg$n_imputations)], vs$pen_vars),
    tau0 = horseshoe_tau0(rc$ref, length(vs$model_pen_vars), cfg),
    gh = gauss_hermite(cfg$gh_nodes),
    grid = time_grid(cfg)
  ))
}
