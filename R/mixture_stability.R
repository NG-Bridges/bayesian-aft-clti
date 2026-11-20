#' Posterior-mean class allocation probabilities (n x 2), rebuilt from stored relabelled mixture draws.
#' The gate is sum-to-zero with K = 2, so gamma[, 1] = gate_logOR / 2 and gamma[, 2] = -gate_logOR / 2.
allocation_probs <- function(m, d, ecol, ms) {
  k_n <- ms$cfg$K_mix
  xp <- design_mat(d, ms$model_pen_vars)
  xd <- design_mat(d, ms$design_vars)
  xw <- design_mat(d, ms$gate_vars)
  logt <- log(d$cr_time)
  ev <- as.integer(d[[ecol]])
  gcol <- sprintf("gate_logOR_%s", ms$gate_vars)
  if (!all(gcol %in% colnames(m))) stop("gate columns missing from stored mixture draws")
  z <- matrix(0, nrow(d), k_n)
  for (s in seq_len(nrow(m))) {
    half <- m[s, gcol] / 2
    gate <- cbind(xw %*% half, xw %*% (-half))
    gate <- gate - apply(gate, 1, max)
    logpi <- gate - log(rowSums(exp(gate)))
    lpc <- matrix(0, nrow(d), k_n)
    for (k in seq_len(k_n)) {
      lp <- as.vector(m[s, sprintf("alpha_c%d", k)] + xp %*% m[s, sprintf("%s_c%d", ms$model_pen_vars, k)] +
        xd %*% m[s, sprintf("%s_c%d", ms$design_vars, k)])
      kk <- m[s, sprintf("k_shape_c%d", k)]
      e <- kk * (logt - lp)
      lpc[, k] <- logpi[, k] - exp(e) + ev * (log(kk) + e - logt)
    }
    ee <- exp(lpc - apply(lpc, 1, max))
    z <- z + ee / rowSums(ee)
  }
  z / nrow(m)
}

#' Adjusted Rand index; NA when either partition puts everyone in one class.
adjusted_rand <- function(a, b) {
  t <- table(a, b)
  if (nrow(t) < 2 || ncol(t) < 2) {
    return(NA_real_)
  }
  s <- sum(choose(t, 2))
  sa <- sum(choose(rowSums(t), 2))
  sb <- sum(choose(colSums(t), 2))
  e <- sa * sb / choose(sum(t), 2)
  den <- (sa + sb) / 2 - e
  if (!is.finite(den) || abs(den) < 1e-12) {
    return(NA_real_)
  }
  (s - e) / den
}

#' Reproducibility of the mixture's hard partition across imputed datasets, for one cause.
#' High ARI means the partition is reproducible, not that the classes are real; read it with the
#' crispness diagnostics (mean_max_prob, mean_entropy) in phase1_fit_diagnostics.csv.
#' @return List: summary (one row), pairs (pairwise ARI), stability (per-patient share in class 1).
mixture_stability <- function(cause, ckpt_dir, diag, ms) {
  cfg <- ms$cfg
  ids <- seq_len(cfg$n_imputations)
  n <- nrow(ms$base)
  ecol <- cause_codes()[[cause]]$event
  lab <- matrix(NA_integer_, n, length(ids))
  pmax_store <- matrix(NA_real_, n, length(ids))
  crisp <- numeric(length(ids))
  for (j in seq_along(ids)) {
    z <- allocation_probs(
      phase1_draws(ckpt_dir, "AFT_mixture", cause, ids[j]),
      prepare_design(ms$imps[[ids[j]]], ms, ms$fix_mom, cfg), ecol, ms
    )
    lab[, j] <- max.col(z)
    pmax_store[, j] <- apply(z, 1, max)
    crisp[j] <- mean(pmax_store[, j] > cfg$ari_prob_cut)
  }
  dg <- diag[diag$model == "AFT_mixture" & diag$cause == cause, ]
  mi <- match(ids, dg$imputation)
  keep <- !is.na(mi) & dg$rhat_lp[mi] < cfg$ari_rhat_lp_max

  pair_ari <- function(cols) {
    if (length(cols) < 2) {
      return(NA_real_)
    }
    apply(combn(cols, 2), 2, function(p) adjusted_rand(lab[, p[1]], lab[, p[2]]))
  }
  a <- pair_ari(seq_along(ids))
  ak <- pair_ari(which(keep))
  q <- quantile(a, c(0, .25, .5, .75, 1), na.rm = TRUE)
  # Agreement with the majority partition: how far one dataset departs from what the 50 agree on
  cons <- ifelse(rowMeans(lab == 1) >= 0.5, 1L, 2L)
  ac <- vapply(seq_along(ids), function(j) adjusted_rand(lab[, j], cons), numeric(1))
  pstab <- rowMeans(lab == 1)
  list(
    summary = data.frame(
      cause = cause, n_datasets = length(ids), n_converged = sum(keep), n_pairs = sum(is.finite(a)),
      ari_median = q[3], ari_q25 = q[2], ari_q75 = q[4], ari_min = q[1], ari_max = q[5],
      ari_median_converged = if (sum(is.finite(ak))) median(ak, na.rm = TRUE) else NA_real_,
      ari_vs_consensus_median = median(ac, na.rm = TRUE),
      frac_prob_gt_cut = mean(crisp),
      median_assigned_prob = median(pmax_store),
      frac_stable_patients = mean(pstab > 0.9 | pstab < 0.1)
    ),
    pairs = data.frame(pair_ari = a),
    stability = data.frame(row = seq_len(n), prop_class1 = pstab)
  )
}
