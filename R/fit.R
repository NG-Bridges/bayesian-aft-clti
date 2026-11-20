.stan_cache <- new.env(parent = emptyenv())

#' Compile a model in `models/` once per session.
stan_model_for <- function(model) {
  name <- c(AFT_single = "aft_single", AFT_frailty = "aft_frailty", AFT_mixture = "aft_mixture")[[model]]
  if (is.null(.stan_cache[[name]])) {
    .stan_cache[[name]] <- rstan::stan_model(file = here::here("models", paste0(name, ".stan")), model_name = name)
  }
  .stan_cache[[name]]
}

stan_data <- function(model, df, ecol, cause, ms) {
  cfg <- ms$cfg
  sd0 <- list(
    N = nrow(df), P = length(ms$model_pen_vars), Q = length(ms$design_vars),
    time = df$cr_time, status = as.integer(df[[ecol]]),
    X = design_mat(df, ms$model_pen_vars), D = design_mat(df, ms$design_vars),
    intercept_scale = cfg$intercept_scale, shape_scale = cfg$shape_scale,
    design_scale = cfg$design_scale, tau_0 = ms$tau0[[cause]],
    slab_scale = cfg$slab_scale, slab_df = cfg$slab_df
  )
  if (model == "AFT_frailty") {
    sd0 <- c(sd0, list(G = cfg$gh_nodes, gh_x = ms$gh$x, gh_logw = ms$gh$logw, frailty_scale = cfg$frailty_scale))
  }
  if (model == "AFT_mixture") {
    sd0 <- c(sd0, list(K = cfg$K_mix, Pw = length(ms$gate_vars), W = design_mat(df, ms$gate_vars), gate_scale = cfg$gate_scale))
  }
  sd0
}

run_stan <- function(model, df, ecol, cause, ms) {
  cfg <- ms$cfg
  rstan::sampling(stan_model_for(model),
    data = stan_data(model, df, ecol, cause, ms),
    chains = cfg$stan_chains, cores = cfg$stan_cores,
    warmup = cfg$stan_warmup, iter = cfg$stan_warmup + cfg$stan_sample,
    seed = cfg$seed, refresh = 0,
    control = list(adapt_delta = cfg$stan_adapt_delta, max_treedepth = cfg$stan_treedepth)
  )
}

converged <- function(g, cfg) {
  g$max_rhat < cfg$rhat_max & g$min_ess > cfg$ess_min & g$div_frac <= cfg$div_frac_max & g$treedepth_hits == 0
}

#' Per-fit convergence on the identified parameters.
#' For the mixture these are on unrelabelled draws; Phase 1 replaces them with relabelled values.
fit_diagnostics <- function(fit, model, cfg) {
  pars <- switch(model,
    AFT_single = c("alpha", "beta", "beta_des", "k_shape", "tau"),
    AFT_frailty = c("alpha", "beta", "beta_des", "k_shape", "tau", "sigma_u"),
    AFT_mixture = c("alpha", "beta", "beta_des", "k_shape", "tau", "gamma_full")
  )
  ss <- rstan::summary(fit, pars = pars)$summary
  lp <- rstan::summary(fit, pars = "lp__")$summary
  ndiv <- rstan::get_num_divergent(fit)
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  td <- sum(vapply(sp, function(x) sum(x[, "treedepth__"] >= cfg$stan_treedepth), numeric(1)))
  g <- data.frame(
    max_rhat = max(ss[, "Rhat"], na.rm = TRUE), min_ess = min(ss[, "n_eff"], na.rm = TRUE),
    rhat_lp = lp[1, "Rhat"], ess_lp = lp[1, "n_eff"],
    n_divergent = ndiv, div_frac = ndiv / (cfg$stan_chains * cfg$stan_sample),
    treedepth_hits = td
  )
  g$converged <- converged(g, cfg)
  g
}

thin_index <- function(nd, n_keep) if (nd > n_keep) round(seq(1, nd, length.out = n_keep)) else seq_len(nd)

#' Reportable draws (single-component and frailty models) as one named matrix.
reportable_draws <- function(fit, model, ms) {
  po <- rstan::extract(fit)
  keep <- thin_index(length(po$alpha), ms$cfg$store_draws)
  m <- cbind(alpha = po$alpha[keep], k_shape = po$k_shape[keep], tau = po$tau[keep])
  if (model == "AFT_frailty") m <- cbind(m, sigma_u = po$sigma_u[keep])
  b <- po$beta[keep, , drop = FALSE]
  colnames(b) <- ms$model_pen_vars
  bd <- po$beta_des[keep, , drop = FALSE]
  colnames(bd) <- ms$design_vars
  cbind(m, b, bd)
}

#' The compact draws kept for prediction (no log_lik, no N-vectors).
prediction_draws <- function(fit, model, n_use) {
  po <- rstan::extract(fit)
  mix <- model == "AFT_mixture"
  keep <- thin_index(if (mix) dim(po$alpha)[1] else length(po$alpha), n_use)
  out <- list(model = model, ndraw = length(keep))
  if (mix) {
    out$alpha <- po$alpha[keep, , drop = FALSE]
    out$k_shape <- po$k_shape[keep, , drop = FALSE]
    out$beta <- po$beta[keep, , , drop = FALSE]
    out$beta_des <- po$beta_des[keep, , , drop = FALSE]
    out$gamma <- po$gamma_full[keep, , , drop = FALSE]
  } else {
    out$alpha <- po$alpha[keep]
    out$k_shape <- po$k_shape[keep]
    out$beta <- po$beta[keep, , drop = FALSE]
    out$beta_des <- po$beta_des[keep, , drop = FALSE]
    if (model == "AFT_frailty") out$sigma_u <- po$sigma_u[keep]
  }
  out
}

#' Relabel mixture draws (Stephens within a fit, then an eGFR anchor across fits) and summarise them.
#' Allocation probabilities are rebuilt in R from thinned draws, so no draws x N x K array leaves Stan.
mixture_postprocess <- function(fit, d, ecol, ms) {
  cfg <- ms$cfg
  K <- cfg$K_mix
  pn <- length(ms$model_pen_vars)
  qn <- length(ms$design_vars)
  pw <- length(ms$gate_vars)
  a <- as.array(fit, pars = c("alpha", "k_shape", "beta", "beta_des", "gamma_full"))
  s0 <- dim(a)[1]
  n_chain <- dim(a)[2]
  a <- a[seq(1, s0, by = max(1, floor(s0 * n_chain / cfg$relabel_draws))), , , drop = FALSE]
  s_n <- dim(a)[1]
  nd <- s_n * n_chain
  nm <- dimnames(a)[[3]]
  # draw index d = s + (c - 1) * S, so array(M, c(S, C, .)) reverses this exactly
  m <- matrix(a, nrow = nd)
  colnames(m) <- nm
  rm(a)

  ia <- match(sprintf("alpha[%d]", seq_len(K)), nm)
  ik <- match(sprintf("k_shape[%d]", seq_len(K)), nm)
  ib <- matrix(match(sprintf("beta[%d,%d]", rep(seq_len(K), pn), rep(seq_len(pn), each = K)), nm), K, pn)
  id <- matrix(match(sprintf("beta_des[%d,%d]", rep(seq_len(K), qn), rep(seq_len(qn), each = K)), nm), K, qn)
  ig <- matrix(match(sprintf("gamma_full[%d,%d]", rep(seq_len(pw), K), rep(seq_len(K), each = pw)), nm), pw, K)
  if (anyNA(c(ia, ik, ib, id, ig))) stop("Could not map mixture parameter names from as.array().")

  xp <- design_mat(d, ms$model_pen_vars)
  xd <- design_mat(d, ms$design_vars)
  xw <- design_mat(d, ms$gate_vars)
  logt <- log(d$cr_time)
  ev <- as.integer(d[[ecol]])
  n <- nrow(d)
  z <- array(0, c(nd, n, K))
  for (s in seq_len(nd)) {
    gate <- xw %*% matrix(m[s, ig], pw, K)
    gate <- gate - apply(gate, 1, max)
    logpi <- gate - log(rowSums(exp(gate)))
    lpc <- matrix(0, n, K)
    for (k in seq_len(K)) {
      lp <- as.vector(m[s, ia[k]] + xp %*% m[s, ib[k, ]] + xd %*% m[s, id[k, ]])
      kk <- m[s, ik[k]]
      e <- kk * (logt - lp)
      lpc[, k] <- logpi[, k] - exp(e) + ev * (log(kk) + e - logt)
    }
    ee <- exp(lpc - apply(lpc, 1, max))
    z[s, , ] <- ee / rowSums(ee)
  }

  st <- try(label.switching::stephens(z), silent = TRUE)
  relabel_ok <- !inherits(st, "try-error") && !is.null(st$permutations)
  perms <- if (relabel_ok) st$permutations else matrix(rep(seq_len(K), each = nd), nd, K)
  ma <- m
  for (s in seq_len(nd)) {
    pp <- perms[s, ]
    if (all(pp == seq_len(K))) next
    ma[s, ia] <- m[s, ia[pp]]
    ma[s, ik] <- m[s, ik[pp]]
    for (k in seq_len(K)) {
      ma[s, ib[k, ]] <- m[s, ib[pp[k], ]]
      ma[s, id[k, ]] <- m[s, id[pp[k], ]]
    }
    ma[s, ig] <- m[s, as.vector(ig[, pp])]
    z[s, , ] <- z[s, , pp]
  }
  switch_frac <- mean(apply(perms, 1, function(pp) any(pp != seq_len(K))))
  rm(m)

  zbar <- apply(z, c(2, 3), mean)
  rm(z)
  if (cfg$anchor_var %in% ms$model_pen_vars) {
    av <- xp[, cfg$anchor_var]
    cls <- vapply(seq_len(K), function(k) sum(zbar[, k] * av) / sum(zbar[, k]), numeric(1))
    anchor_used <- cfg$anchor_var
  } else {
    cls <- colMeans(ma[, ia, drop = FALSE])
    anchor_used <- "alpha"
  }
  ap <- order(cls)
  mb <- ma
  if (!all(ap == seq_len(K))) {
    mb[, ia] <- ma[, ia[ap]]
    mb[, ik] <- ma[, ik[ap]]
    for (k in seq_len(K)) {
      mb[, ib[k, ]] <- ma[, ib[ap[k], ]]
      mb[, id[k, ]] <- ma[, id[ap[k], ]]
    }
    mb[, ig] <- ma[, as.vector(ig[, ap])]
    zbar <- zbar[, ap, drop = FALSE]
  }
  rm(ma)

  # Raw R-hat is inflated by label switching, so convergence is judged on the relabelled draws
  arr <- array(mb, c(s_n, n_chain, ncol(mb)), dimnames = list(NULL, NULL, colnames(mb)))
  mon <- rstan::monitor(arr, warmup = 0, print = FALSE)
  rm(arr)

  keep <- thin_index(nd, cfg$store_draws)
  cols <- list()
  for (k in seq_len(K)) {
    cols[[sprintf("alpha_c%d", k)]] <- mb[keep, ia[k]]
    cols[[sprintf("k_shape_c%d", k)]] <- mb[keep, ik[k]]
    for (j in seq_len(pn)) cols[[sprintf("%s_c%d", ms$model_pen_vars[j], k)]] <- mb[keep, ib[k, j]]
    for (j in seq_len(qn)) cols[[sprintf("%s_c%d", ms$design_vars[j], k)]] <- mb[keep, id[k, j]]
  }
  cols[["alpha_c2_minus_c1"]] <- mb[keep, ia[2]] - mb[keep, ia[1]]
  for (v in seq_len(pw)) cols[[sprintf("gate_logOR_%s", ms$gate_vars[v])]] <- mb[keep, ig[v, 1]] - mb[keep, ig[v, 2]]

  ent <- -rowSums(zbar * log(pmax(zbar, 1e-12))) / log(K)
  list(
    draws = do.call(cbind, cols),
    diag = data.frame(
      max_rhat_rel = max(mon$Rhat, na.rm = TRUE),
      min_ess_rel = min(mon$Bulk_ESS, na.rm = TRUE),
      relabel_ok = relabel_ok,
      anchor_used = anchor_used,
      switch_frac = switch_frac,
      mean_max_prob = mean(apply(zbar, 1, max)),
      mean_entropy = mean(ent),
      class1_share = mean(zbar[, 1]),
      stringsAsFactors = FALSE
    )
  )
}

.diag_cols <- c(
  "max_rhat", "min_ess", "rhat_lp", "ess_lp", "n_divergent", "div_frac",
  "treedepth_hits", "converged", "cor_sigmau_kshape", "sigma_u_mean",
  "max_rhat_raw", "min_ess_raw", "relabel_ok", "anchor_used", "switch_frac",
  "mean_max_prob", "mean_entropy", "class1_share"
)

#' Pad diagnostics to the shared columns so all three models rbind.
pad_diag <- function(dg) {
  for (cl in setdiff(.diag_cols, names(dg))) dg[[cl]] <- if (cl == "anchor_used") NA_character_ else NA_real_
  dg[, .diag_cols, drop = FALSE]
}

#' One Phase 1 fit reduced to what is kept: information criteria, diagnostics, reportable draws.
summarise_inference_fit <- function(fit, model, d, ecol, ms) {
  ll <- rstan::extract(fit, pars = "log_lik")$log_lik
  ic <- ic_from_loglik(ll)
  rm(ll)
  dg <- fit_diagnostics(fit, model, ms$cfg)
  if (model == "AFT_mixture") {
    mp <- mixture_postprocess(fit, d, ecol, ms)
    dg$max_rhat_raw <- dg$max_rhat
    dg$min_ess_raw <- dg$min_ess
    dg$max_rhat <- mp$diag$max_rhat_rel
    dg$min_ess <- mp$diag$min_ess_rel
    dg$converged <- converged(dg, ms$cfg)
    dg <- cbind(dg, mp$diag[, setdiff(names(mp$diag), c("max_rhat_rel", "min_ess_rel"))])
    drw <- mp$draws
  } else {
    if (model == "AFT_frailty") {
      po <- rstan::extract(fit, pars = c("sigma_u", "k_shape"))
      dg$cor_sigmau_kshape <- cor(po$sigma_u, po$k_shape)
      dg$sigma_u_mean <- mean(po$sigma_u)
    }
    drw <- reportable_draws(fit, model, ms)
  }
  list(ic = ic, diag = pad_diag(dg), draws = drw)
}
