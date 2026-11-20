#' Out-of-fold CIFs for one model and dataset: saved Phase 3 predictions (AFT) or the RSF checkpoint.
#' @return List: amputation, death (n x horizons), or NULL if not available.
load_oof <- function(model, i, ckpt_dir, cfg) {
  if (model == "RSF") {
    f <- checkpoint_path(ckpt_dir, "rsf_cv", i)
    if (!file.exists(f)) {
      return(NULL)
    }
    r <- readRDS(f)
    return(list(amputation = as.matrix(r$amputation), death = as.matrix(r$death)))
  }
  f <- results_path("cv_predictions", basename(ckpt_dir), sprintf("oof_%s_grid%d_%d.rds", model, cfg$grid_n, i))
  if (!file.exists(f)) {
    return(NULL)
  }
  readRDS(f)
}

#' Aalen-Johansen CIF of both causes at `times`.
#' @return Matrix: length(times) x 2 (amputation, death).
aj_cif <- function(time, status, times) {
  o <- order(time)
  time <- time[o]
  status <- status[o]
  ut <- sort(unique(time[status > 0]))
  nrisk <- length(time) - findInterval(ut, time, left.open = TRUE)
  d1 <- tabulate(match(time[status == 1], ut), length(ut))
  d2 <- tabulate(match(time[status == 2], ut), length(ut))
  s_minus <- c(1, cumprod(1 - (d1 + d2) / nrisk))[seq_along(ut)]
  j <- findInterval(times, ut)
  c1 <- c(0, cumsum(s_minus * d1 / nrisk))
  c2 <- c(0, cumsum(s_minus * d2 / nrisk))
  cbind(amputation = c1[j + 1], death = c2[j + 1])
}

#' Jackknife pseudo-values of the Aalen-Johansen CIF (Andersen and Perme 2010), checked against prodlim.
pseudo_values <- function(base, cfg) {
  tt <- base$cr_time
  st <- base$cr_status
  n <- length(tt)
  full <- aj_cif(tt, st, cfg$horizons)
  pov <- list(amputation = matrix(NA_real_, n, length(cfg$horizons)), death = matrix(NA_real_, n, length(cfg$horizons)))
  for (j in seq_len(n)) {
    loo <- aj_cif(tt[-j], st[-j], cfg$horizons)
    pov$amputation[j, ] <- n * full[, 1] - (n - 1) * loo[, 1]
    pov$death[j, ] <- n * full[, 2] - (n - 1) * loo[, 2]
  }
  pl <- prodlim::prodlim(Hist(cr_time, cr_status) ~ 1, data = base)
  for (cs in names(cause_codes())) {
    ref <- as.numeric(unlist(predict(pl, times = cfg$horizons, cause = cause_codes()[[cs]]$code)))
    if (max(abs(ref - full[, cs])) > 1e-8) stop("Aalen-Johansen estimate disagrees with prodlim for ", cs)
  }
  pov
}

cll <- function(p) log(-log(1 - p))
clamp_p <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
smooth_fit <- function(p, pov, df) lm(pov ~ splines::ns(cll(p), df = df))

#' Calibration of predicted risk `p` against pseudo-values: O/E, cloglog intercept and slope, ICI/E50/E90.
#' The quasi family with constant variance is used because pseudo-values can leave [0, 1].
calib_metrics <- function(p, pov, df) {
  p <- clamp_p(p)
  fam <- quasi(link = "cloglog", variance = "constant")
  safe <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) NULL)
  err <- abs(fitted(smooth_fit(p, pov, df)) - p)
  g0 <- safe(glm(pov ~ 1 + offset(cll(p)), family = fam, start = 0, control = glm.control(maxit = 100)))
  g1 <- safe(glm(pov ~ cll(p), family = fam, start = c(0, 1), control = glm.control(maxit = 100)))
  c(
    OE = mean(pov) / mean(p),
    intercept = if (!is.null(g0) && g0$converged) unname(coef(g0)[1]) else NA_real_,
    slope = if (!is.null(g1) && g1$converged) unname(coef(g1)[2]) else NA_real_,
    ICI = mean(err), E50 = unname(median(err)), E90 = unname(quantile(err, 0.9)),
    mean_pred = mean(p), observed = mean(pov)
  )
}

#' Rebuilt predictions must reproduce the Phase 2/3 scores; otherwise they are not the reported model.
verify_oof <- function(oof, avail, imps, cfg) {
  saved <- list(
    AFT_single = results_path("phase3_bayes_cv.rds"),
    RSF = results_path("phase2_rsf_cv.rds")
  )
  rows <- list()
  for (mdl in cfg$calib_models) {
    if (!file.exists(saved[[mdl]])) stop("Missing ", saved[[mdl]], "; run scripts/02_fit_models.R first.")
    sv <- readRDS(saved[[mdl]])
    sv <- sv[sv$model == mdl, ]
    for (i in head(avail, cfg$calib_verify_n)) {
      for (cs in names(cause_codes())) {
        d <- imps[[i]]
        rb <- score_cr(d$cr_time, d$cr_status, cause_codes()[[cs]]$code, oof[[mdl]][[i]][[cs]], cfg$horizons)
        z <- sv[sv$cause == cs & sv$imputation == i, ]
        rows[[length(rows) + 1]] <- data.frame(
          model = mdl, cause = cs, imputation = i, months = cfg$horizons_months,
          AUC_rebuilt = rb$AUC, AUC_saved = z$AUC[match(cfg$horizons_months, z$months)],
          Brier_rebuilt = rb$Brier, Brier_saved = z$Brier[match(cfg$horizons_months, z$months)]
        )
      }
    }
  }
  ver <- do.call(rbind, rows)
  ver$dAUC <- ver$AUC_rebuilt - ver$AUC_saved
  ver$dBrier <- ver$Brier_rebuilt - ver$Brier_saved
  if (any(abs(ver$dAUC) > 1e-6 | abs(ver$dBrier) > 1e-6, na.rm = TRUE)) {
    stop("Out-of-fold predictions do not reproduce the saved Phase 2/3 scores.")
  }
  ver
}

#' Calibration per dataset (imputation spread) and on predictions averaged over datasets (primary, bootstrap CIs).
#' The bootstrap is conditional on the single fold partition shared by every dataset.
calibration_analysis <- function(oof, pov, avail, cfg) {
  h_n <- length(cfg$horizons)
  metric_names <- c("OE", "intercept", "slope", "ICI", "E50", "E90", "mean_pred", "observed")
  causes <- names(cause_codes())

  bd_rows <- list()
  for (mdl in cfg$calib_models) {
    for (i in avail) {
      for (cs in causes) {
        for (h in seq_len(h_n)) {
          m <- calib_metrics(oof[[mdl]][[i]][[cs]][, h], pov[[cs]][, h], cfg$calib_spline_df)
          bd_rows[[length(bd_rows) + 1]] <- data.frame(model = mdl, cause = cs, imputation = i, months = cfg$horizons_months[h], t(m))
        }
      }
    }
  }
  by_ds <- do.call(rbind, bd_rows)
  spread <- do.call(rbind, lapply(split(by_ds, list(by_ds$model, by_ds$cause, by_ds$months), drop = TRUE), function(z) {
    out <- data.frame(model = z$model[1], cause = z$cause[1], months = z$months[1])
    for (mn in metric_names) {
      out[[paste0(mn, "_mean")]] <- mean(z[[mn]], na.rm = TRUE)
      out[[paste0(mn, "_p025")]] <- unname(quantile(z[[mn]], 0.025, na.rm = TRUE))
      out[[paste0(mn, "_p975")]] <- unname(quantile(z[[mn]], 0.975, na.rm = TRUE))
    }
    out
  }))
  spread <- spread[order(spread$cause, spread$months, spread$model), ]

  pbar <- lapply(oof, function(l) {
    lapply(setNames(causes, causes), function(cs) Reduce(`+`, lapply(l[avail], `[[`, cs)) / length(avail))
  })
  n <- nrow(pov[[1]])
  set.seed(cfg$seed)
  boot_idx <- replicate(cfg$calib_n_boot, sample.int(n, n, replace = TRUE), simplify = FALSE)

  sum_rows <- list()
  dec_rows <- list()
  crv_rows <- list()
  for (mdl in cfg$calib_models) {
    for (cs in causes) {
      for (h in seq_len(h_n)) {
        p <- clamp_p(pbar[[mdl]][[cs]][, h])
        pv <- pov[[cs]][, h]
        est <- calib_metrics(p, pv, cfg$calib_spline_df)
        pg <- seq(quantile(p, 0.01), quantile(p, 0.99), length.out = 100)
        nd <- data.frame(p = pg)
        crv <- predict(smooth_fit(p, pv, cfg$calib_spline_df), newdata = nd)
        b_m <- matrix(NA_real_, cfg$calib_n_boot, length(metric_names), dimnames = list(NULL, metric_names))
        b_c <- matrix(NA_real_, cfg$calib_n_boot, length(pg))
        for (b in seq_len(cfg$calib_n_boot)) {
          ix <- boot_idx[[b]]
          b_m[b, ] <- calib_metrics(p[ix], pv[ix], cfg$calib_spline_df)
          # Each bootstrap smoother is evaluated only inside its own data range
          inr <- pg >= min(p[ix]) & pg <= max(p[ix])
          b_c[b, inr] <- predict(smooth_fit(p[ix], pv[ix], cfg$calib_spline_df), newdata = nd[inr, , drop = FALSE])
        }
        row <- data.frame(model = mdl, cause = cs, months = cfg$horizons_months[h])
        for (mn in metric_names) {
          row[[mn]] <- unname(est[mn])
          row[[paste0(mn, "_lo")]] <- unname(quantile(b_m[, mn], 0.025, na.rm = TRUE))
          row[[paste0(mn, "_hi")]] <- unname(quantile(b_m[, mn], 0.975, na.rm = TRUE))
        }
        sum_rows[[length(sum_rows) + 1]] <- row
        crv_rows[[length(crv_rows) + 1]] <- data.frame(
          model = mdl, cause = cs, months = cfg$horizons_months[h], predicted = pg, observed = crv,
          lo = apply(b_c, 2, quantile, 0.025, na.rm = TRUE), hi = apply(b_c, 2, quantile, 0.975, na.rm = TRUE)
        )
        grp <- cut(rank(p, ties.method = "first"), cfg$calib_n_groups, labels = FALSE)
        dec_rows[[length(dec_rows) + 1]] <- do.call(rbind, lapply(seq_len(cfg$calib_n_groups), function(g) {
          ii <- grp == g
          data.frame(
            model = mdl, cause = cs, months = cfg$horizons_months[h], group = g, n = sum(ii),
            predicted = mean(p[ii]), observed = mean(pv[ii]), se = sd(pv[ii]) / sqrt(sum(ii))
          )
        }))
      }
    }
  }
  summ <- do.call(rbind, sum_rows)
  dec <- do.call(rbind, dec_rows)
  dec$lo <- dec$observed - 1.96 * dec$se
  dec$hi <- dec$observed + 1.96 * dec$se
  summ <- summ[order(summ$cause, summ$months, summ$model), ]
  list(by_dataset = by_ds, spread = spread, summary = summ, deciles = dec, curves = do.call(rbind, crv_rows), pbar = pbar)
}
