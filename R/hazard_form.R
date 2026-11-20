#' Rubin-pooled multivariate Wald test (D1; van Buuren, Flexible Imputation of Missing Data, 5.3.2)
#' of H0: theta = theta0 for the parameters in `ests` (one vector per imputation).
pool_d1 <- function(ests, vcovs, theta0 = NULL) {
  m <- length(ests)
  k <- length(ests[[1]])
  if (is.null(theta0)) theta0 <- rep(0, k)
  qbar <- Reduce(`+`, ests) / m
  ubar <- Reduce(`+`, vcovs) / m
  dev <- lapply(ests, function(q) matrix(q - qbar, ncol = 1))
  b <- Reduce(`+`, lapply(dev, function(v) v %*% t(v))) / (m - 1)
  uinv <- tryCatch(solve(ubar), error = function(e) NULL)
  if (is.null(uinv)) {
    return(NULL)
  }
  r1 <- (1 + 1 / m) * sum(diag(b %*% uinv)) / k
  d <- matrix(qbar - theta0, ncol = 1)
  d1 <- as.numeric(t(d) %*% uinv %*% d) / (k * (1 + r1))
  t_ <- k * (m - 1)
  v1 <- if (t_ > 4) 4 + (t_ - 4) * (1 + (1 - 2 / t_) / r1)^2 else t_ * (1 + 1 / k) * (1 + 1 / r1)^2 / 2
  list(
    stat = d1, df1 = k, df2 = v1, p = pf(d1, k, v1, lower.tail = FALSE),
    r1 = r1, fmi = (r1 + 2 / (v1 + 3)) / (r1 + 1),
    est = qbar, se = sqrt(diag(ubar + (1 + 1 / m) * b))
  )
}

#' Cox-Snell residual slope of a survreg Weibull fit: 1 means hazard accumulates at the right rate.
cox_snell_slope <- function(fit, d, ecol) {
  lp <- predict(fit, newdata = d, type = "lp")
  cs <- (d$cr_time / exp(lp))^(1 / fit$scale)
  sf <- survival::survfit(survival::Surv(cs, d[[ecol]]) ~ 1)
  h <- -log(pmax(sf$surv, 1e-12))
  x <- sf$time
  ok <- is.finite(h) & x > 0
  if (sum(ok) < 5) {
    return(NA_real_)
  }
  unname(coef(lm(h[ok] ~ 0 + x[ok])))
}

#' Weibull adequacy per cause on the raw completed datasets: Cox-Snell slope, and a generalized gamma
#' likelihood-ratio test with the shape Q pooled across imputations (Q = 1 is the Weibull).
#' @param core `P$core` from the prepared data; the predictor order matches the original analysis.
#' @return List: pooled (one row per test and cause), per_draw (per-imputation LRT).
hazard_form_tests <- function(imps, core, cfg) {
  design <- c("ARM", "FACTOR_B", if (cfg$include_cohort) "COHORT")
  vars <- c(setdiff(core, c("ARM", "FACTOR_B")), "GNRI", design)
  form <- as.formula(paste("Surv(cr_time, .ev) ~", paste(vars, collapse = " + ")))
  pooled <- list()
  per_draw <- list()
  for (cause in names(cause_codes())) {
    slopes <- c()
    qs <- list()
    qv <- list()
    gg <- list()
    for (i in seq_len(cfg$n_imputations)) {
      d <- imps[[i]]
      d$.ev <- as.integer(d$cr_status == cause_codes()[[cause]]$code)
      fw <- tryCatch(survival::survreg(form, data = d, dist = "weibull"), error = function(e) NULL)
      if (is.null(fw)) next
      slopes <- c(slopes, cox_snell_slope(fw, d, ".ev"))
      fg <- tryCatch(flexsurv::flexsurvreg(form, data = d, dist = "gengamma"), error = function(e) NULL)
      fw2 <- tryCatch(flexsurv::flexsurvreg(form, data = d, dist = "weibull"), error = function(e) NULL)
      if (is.null(fg) || is.null(fw2)) next
      st <- max(2 * (fg$loglik - fw2$loglik), 0)
      gg[[length(gg) + 1]] <- data.frame(
        draw = i, family = "hazard_form", test = "weibull_vs_gengamma", cause = cause,
        chisq = st, df = 1, p = pchisq(st, 1, lower.tail = FALSE), stringsAsFactors = FALSE
      )
      qe <- fg$res.t["Q", "est"]
      qse <- fg$res.t["Q", "se"]
      if (is.finite(qe) && is.finite(qse)) {
        qs[[length(qs) + 1]] <- c(Q = qe)
        qv[[length(qv) + 1]] <- matrix(qse^2, 1, 1, dimnames = list("Q", "Q"))
      }
    }
    pooled[[length(pooled) + 1]] <- data.frame(
      family = "hazard_form", test = "cox_snell_slope_weibull", cause = cause, k_df = NA,
      n_draws = sum(is.finite(slopes)), p_min_draw = NA, p_max_draw = NA, n_sig05 = NA,
      D1 = mean(slopes, na.rm = TRUE), df1 = NA, df2 = NA, p_pooled = NA,
      riv = NA, fmi = NA, est1 = mean(slopes, na.rm = TRUE),
      se1 = sd(slopes, na.rm = TRUE), mde_TR = NA, aic_small = NA, aic_large = NA,
      stringsAsFactors = FALSE
    )
    per_draw <- c(per_draw, gg)
    pd <- if (length(qs) >= 2) pool_d1(qs, qv, theta0 = 1) else NULL
    ps <- vapply(gg, function(r) r$p, numeric(1))
    pooled[[length(pooled) + 1]] <- data.frame(
      family = "hazard_form", test = "weibull_vs_gengamma", cause = cause, k_df = 1,
      n_draws = length(gg), p_min_draw = min(ps), p_max_draw = max(ps), n_sig05 = sum(ps < 0.05),
      D1 = if (is.null(pd)) NA else pd$stat, df1 = if (is.null(pd)) NA else pd$df1,
      df2 = if (is.null(pd)) NA else pd$df2, p_pooled = if (is.null(pd)) NA else pd$p,
      riv = if (is.null(pd)) NA else pd$r1, fmi = if (is.null(pd)) NA else pd$fmi,
      est1 = if (is.null(pd)) NA else pd$est[1], se1 = if (is.null(pd)) NA else pd$se[1],
      # Smallest time ratio detectable at 80% power given the pooled SE
      mde_TR = if (is.null(pd)) NA else exp(2.802 * pd$se[1]),
      aic_small = NA, aic_large = NA, stringsAsFactors = FALSE
    )
  }
  po <- do.call(rbind, pooled)
  po$p_BH <- NA_real_
  ix <- which(is.finite(po$p_pooled))
  po$p_BH[ix] <- p.adjust(po$p_pooled[ix], "BH")
  list(pooled = po, per_draw = do.call(rbind, per_draw))
}
