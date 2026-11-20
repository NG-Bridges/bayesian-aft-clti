#' WAIC and DIC3 (Celeux et al. 2006) from one pointwise log-likelihood matrix (draws x N).
#' DIC3 replaces the plug-in deviance with the posterior-mean predictive density, so it is label-invariant.
ic_from_loglik <- function(ll) {
  n <- ncol(ll)
  cmax <- apply(ll, 2, max)
  lppd_i <- cmax + log(colMeans(exp(sweep(ll, 2, cmax, "-"))))
  pw_i <- apply(ll, 2, var)
  elpd_i <- lppd_i - pw_i
  lppd <- sum(lppd_i)
  dbar <- -2 * mean(rowSums(ll))
  p_dic3 <- dbar + 2 * lppd
  list(
    waic = -2 * sum(elpd_i),
    se_waic = 2 * sqrt(n) * sd(elpd_i),
    p_waic = sum(pw_i),
    lppd = lppd,
    dic3 = dbar + p_dic3,
    p_dic3 = p_dic3,
    Dbar = dbar,
    max_p_waic_i = max(pw_i),
    n_p_waic_gt_04 = sum(pw_i > 0.4),
    elpd_i = elpd_i
  )
}

#' Rubin's rules with the Barnard-Rubin small-sample df; used here only to report FMI.
rubin_pool <- function(ests, vars_diag, nu_com) {
  m <- nrow(ests)
  qbar <- colMeans(ests)
  ubar <- colMeans(vars_diag)
  b <- apply(ests, 2, var)
  tvar <- ubar + (1 + 1 / m) * b
  riv <- ((1 + 1 / m) * b) / ubar
  lambda <- ((1 + 1 / m) * b) / tvar
  lambda[!is.finite(lambda) | lambda < 1e-8] <- 1e-8
  df_old <- (m - 1) / lambda^2
  dfobs <- ((nu_com + 1) / (nu_com + 3)) * nu_com * (1 - lambda)
  df <- df_old * dfobs / (df_old + dfobs)
  se <- sqrt(tvar)
  tval <- qbar / se
  tcrit <- qt(0.975, df)
  data.frame(
    variable = colnames(ests), est = qbar, se = se,
    lo = qbar - tcrit * se, hi = qbar + tcrit * se,
    t = tval, df = df, p = 2 * pt(-abs(tval), df),
    riv = riv, fmi = (riv + 2 / (df + 3)) / (riv + 1), row.names = NULL
  )
}

#' Marginal log-likelihood of the frailty model per patient, recomputed in R to check the quadrature.
gh_loglik <- function(logt, ev, lp, kk, sig, nodes) {
  q <- gauss_hermite(nodes)
  e <- kk * outer(logt - lp, sqrt(2) * sig * q$x, "-")
  lq <- sweep(-exp(e), 2, q$logw, "+") + ev * (log(kk) + e - logt)
  m <- apply(lq, 1, max)
  m + log(rowSums(exp(lq - m)))
}
