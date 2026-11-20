phase1_draws <- function(ckpt_dir, model, cause, i) {
  p <- checkpoint_path(ckpt_dir, "inf", model, cause, i)
  if (!file.exists(p)) stop("Missing Phase 1 checkpoint ", p)
  readRDS(p)$draws
}

summarise_draws <- function(x) {
  c(
    median = unname(median(x)), lo95 = unname(quantile(x, 0.025)), hi95 = unname(quantile(x, 0.975)),
    sd = sd(x), p_direction = max(mean(x > 0), mean(x < 0))
  )
}

#' Stop unless the stacked draws reproduce the pooled Phase 1 summaries (i.e. the same checkpoints).
check_against_pooled <- function(pooled, model, cause, m) {
  z <- pooled[pooled$model == model & pooled$cause == cause & pooled$variable %in% colnames(m), ]
  q <- t(apply(m[, z$variable, drop = FALSE], 2, quantile, c(0.5, 0.025, 0.975), names = FALSE))
  if (max(abs(q - as.matrix(z[, c("post_median", "lo95", "hi95")]))) > 1e-6) {
    stop("Stacked draws do not reproduce phase1_coefficients_pooled.csv for ", model, "/", cause)
  }
}

#' Single-component model, amputation minus death, per coefficient and shape.
#' The two cause-specific fits have separate likelihoods, so pairing draw s of each is a valid joint draw;
#' differencing within dataset before stacking keeps the correlation the shared imputation induces.
between_cause_differences <- function(ckpt_dir, pooled, n_imp) {
  a <- list()
  d <- list()
  dif <- list()
  for (i in seq_len(n_imp)) {
    ai <- phase1_draws(ckpt_dir, "AFT_single", "amputation", i)
    di <- phase1_draws(ckpt_dir, "AFT_single", "death", i)
    v <- intersect(colnames(ai), colnames(di))
    n <- min(nrow(ai), nrow(di))
    a[[i]] <- ai
    d[[i]] <- di
    dif[[i]] <- ai[seq_len(n), v, drop = FALSE] - di[seq_len(n), v, drop = FALSE]
  }
  a <- do.call(rbind, a)
  d <- do.call(rbind, d)
  dif <- do.call(rbind, dif)
  check_against_pooled(pooled, "AFT_single", "amputation", a)
  check_against_pooled(pooled, "AFT_single", "death", d)
  bc <- do.call(rbind, lapply(setdiff(colnames(dif), c("alpha", "tau")), function(v) {
    sa <- summarise_draws(a[, v])
    sd_ <- summarise_draws(d[, v])
    sf <- summarise_draws(dif[, v])
    data.frame(
      variable = v,
      amp_median = sa[["median"]], amp_lo95 = sa[["lo95"]], amp_hi95 = sa[["hi95"]],
      death_median = sd_[["median"]], death_lo95 = sd_[["lo95"]], death_hi95 = sd_[["hi95"]],
      diff_median = sf[["median"]], diff_lo95 = sf[["lo95"]], diff_hi95 = sf[["hi95"]],
      diff_sd = sf[["sd"]], diff_p_direction = sf[["p_direction"]],
      sd_if_independent = sqrt(sa[["sd"]]^2 + sd_[["sd"]]^2),
      sd_ratio = sf[["sd"]] / sqrt(sa[["sd"]]^2 + sd_[["sd"]]^2),
      excludes_zero = sf[["lo95"]] > 0 | sf[["hi95"]] < 0
    )
  }))
  bc[order(-abs(bc$diff_median)), ]
}

#' Gated mixture, Class 2 minus Class 1 within each (relabelled) draw, stacked over datasets.
mixture_class_differences <- function(ckpt_dir, pooled, n_imp) {
  do.call(rbind, lapply(names(cause_codes()), function(cs) {
    m <- do.call(rbind, lapply(seq_len(n_imp), function(i) phase1_draws(ckpt_dir, "AFT_mixture", cs, i)))
    check_against_pooled(pooled, "AFT_mixture", cs, m)
    stems <- sub("_c1$", "", grep("_c1$", colnames(m), value = TRUE))
    stems <- stems[paste0(stems, "_c2") %in% colnames(m) & stems != "alpha"]
    do.call(rbind, lapply(stems, function(v) {
      s <- summarise_draws(m[, paste0(v, "_c2")] - m[, paste0(v, "_c1")])
      data.frame(
        cause = cs, variable = v, median = s[["median"]], lo95 = s[["lo95"]], hi95 = s[["hi95"]],
        sd = s[["sd"]], p_direction = s[["p_direction"]], excludes_zero = s[["lo95"]] > 0 | s[["hi95"]] < 0
      )
    }))
  }))
}
