#' Prediction time grid: quadratic spacing, dense near 0 where the shape < 1 amputation hazard changes fastest.
time_grid <- function(cfg) {
  tg <- sort(unique(c(max(cfg$horizons) * (seq(0, 1, length.out = cfg$grid_n))^2, cfg$horizons)))
  ix <- match(cfg$horizons, tg)
  stopifnot(!anyNA(ix), tg[1] == 0)
  # Column j of `ut` sums the first j interval increments, i.e. gives the CIF at grid point j + 1;
  # `ut_h` keeps only the columns for the horizons, which is all prediction needs
  ut <- outer(seq_len(length(tg) - 1), seq_len(length(tg) - 1), "<=") * 1
  list(tgrid = tg, h_ix = ix, ut = ut, ut_h = ut[, ix - 1, drop = FALSE])
}

#' Cumulative incidence of both causes from a pair of cause-specific survival matrices (n x grid).
#' Each interval's drop in all-cause survival is split by the cumulative-hazard increments,
#' the discrete form of CIF_k = int h_k S, which never evaluates h_k at t = 0.
#' @param ut Summation matrix from `time_grid()`: `ut` for every grid point after 0, `ut_h` for the horizons.
cif_pair <- function(s1, s2, ut) {
  g <- ncol(s1)
  h1 <- -log(pmax(s1, 1e-300))
  h2 <- -log(pmax(s2, 1e-300))
  dh1 <- h1[, -1, drop = FALSE] - h1[, -g, drop = FALSE]
  dh2 <- h2[, -1, drop = FALSE] - h2[, -g, drop = FALSE]
  dh <- dh1 + dh2
  sa <- s1 * s2
  ds <- sa[, -g, drop = FALSE] - sa[, -1, drop = FALSE]
  ds[ds < 0] <- 0
  w1 <- ifelse(dh > 1e-12, dh1 / dh, 0.5)
  list(amp = (ds * w1) %*% ut, death = (ds * (1 - w1)) %*% ut)
}

#' Weibull survival on the grid. S(0) is set to 1 because exp(-k lp) can overflow and Inf * 0 is NaN.
weib_S <- function(lp, kk, tgrid) {
  s <- exp(-(exp(-kk * lp) %o% (tgrid^kk)))
  s[, 1] <- 1
  s
}

#' Cause-specific survival (n x grid) for posterior draw `s` of any of the three models.
surv_draw <- function(dr, s, xp, xd, xw, ms) {
  tgrid <- ms$grid$tgrid
  if (dr$model == "AFT_mixture") {
    K <- ncol(dr$alpha)
    gate <- xw %*% matrix(dr$gamma[s, , ], nrow = ncol(xw), ncol = K)
    gate <- gate - apply(gate, 1, max)
    pim <- exp(gate)
    pim <- pim / rowSums(pim)
    out <- matrix(0, nrow(xp), length(tgrid))
    for (k in seq_len(K)) {
      lp <- as.vector(dr$alpha[s, k] + xp %*% dr$beta[s, k, ] + xd %*% dr$beta_des[s, k, ])
      out <- out + pim[, k] * weib_S(lp, dr$k_shape[s, k], tgrid)
    }
    out[, 1] <- 1
    return(out)
  }
  lp <- as.vector(dr$alpha[s] + xp %*% dr$beta[s, ] + xd %*% dr$beta_des[s, ])
  kk <- dr$k_shape[s]
  if (dr$model == "AFT_frailty") {
    # Population-averaged over the same quadrature as the likelihood;
    # H(t | u) = H0(t) exp(-k u), so each node only rescales H0
    w <- exp(ms$gh$logw)
    h0 <- exp(-kk * lp) %o% (tgrid^kk)
    out <- matrix(0, nrow(xp), length(tgrid))
    for (g in seq_along(ms$gh$x)) out <- out + w[g] * exp(-h0 * exp(-kk * sqrt(2) * dr$sigma_u[s] * ms$gh$x[g]))
    out[, 1] <- 1
    return(out)
  }
  weib_S(lp, kk, tgrid)
}

#' Posterior-mean CIF of both causes at the horizons, pairing draw s of each cause.
predict_cif <- function(dr_amp, dr_dth, xp, xd, xw, ms) {
  nd <- min(dr_amp$ndraw, dr_dth$ndraw)
  a <- matrix(0, nrow(xp), length(ms$grid$h_ix))
  dm <- a
  for (s in seq_len(nd)) {
    cf <- cif_pair(surv_draw(dr_amp, s, xp, xd, xw, ms), surv_draw(dr_dth, s, xp, xd, xw, ms), ms$grid$ut_h)
    a <- a + cf$amp
    dm <- dm + cf$death
  }
  list(amputation = a / nd, death = dm / nd)
}

#' Competing-risks AUC, Brier and IPA of predicted CIFs at each horizon.
#' The formula must use bare Hist(): Score() identifies the response type by matching the call name.
score_cr <- function(time, status, cause, risk, horizons) {
  dd <- data.frame(time = time, status = status)
  sc <- tryCatch(
    riskRegression::Score(list(m = risk),
      formula = Hist(time, status) ~ 1, data = dd,
      cause = cause, times = horizons, metrics = c("AUC", "Brier"), summary = "ipa",
      se.fit = FALSE, conf.int = FALSE, cens.model = "km"
    ),
    error = function(e) {
      message("Score() failed: ", conditionMessage(e))
      NULL
    }
  )
  out <- data.frame(times = horizons, AUC = NA_real_, Brier = NA_real_, IPA = NA_real_)
  if (!is.null(sc)) {
    a <- as.data.frame(sc$AUC$score)
    b <- as.data.frame(sc$Brier$score)
    a <- a[a$model == "m", c("times", "AUC"), drop = FALSE]
    b <- b[b$model == "m", intersect(c("times", "Brier", "IPA"), names(b)), drop = FALSE]
    m <- merge(a, b, by = "times")
    for (nm in c("AUC", "Brier", "IPA")) {
      if (nm %in% names(m)) out[match(m$times, out$times), nm] <- m[[nm]]
    }
  }
  out
}
