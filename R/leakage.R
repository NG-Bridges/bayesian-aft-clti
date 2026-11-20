#' Score predictions within a subset of patients.
score_subset <- function(base, risk, idx, cause, times) {
  score_cr(base$cr_time[idx], base$cr_status[idx], cause, risk[idx, , drop = FALSE], times)
}

#' Does outcome-informed imputation of albumin inflate the cross-validated RSF results?
#' Test A: GNRI's contribution to AUC (full RSF minus an RSF without GNRI), imputed minus measured
#' albumin (difference-in-differences), plus a patient bootstrap on one dataset.
#' Test B: cause-specific Cox with a GNRI x imputed-albumin interaction, Rubin-pooled.
#' @return List of data frames: metrics, did, bootstrap, cox.
leakage_check <- function(ms, fold_id, ckpt_dir, cache_dir) {
  cfg <- ms$cfg
  base <- ms$base
  imps <- ms$imps
  n <- nrow(base)
  hz <- cfg$horizons
  codes <- vapply(cause_codes(), `[[`, integer(1), "code")
  ids <- seq_len(cfg$rsf_n_imputations)

  # base keeps raw albumin, so GNRI is missing exactly where albumin was imputed
  r_imp <- is.na(base$GNRI)
  if (!any(r_imp)) stop("GNRI has no missing values in the base data; imputed patients cannot be identified.")
  groups <- list(all = seq_len(n), observed = which(!r_imp), imputed = which(r_imp))

  oof <- lapply(ids, function(i) readRDS(checkpoint_path(ckpt_dir, "rsf_cv", i)))

  # If GNRI is the only imputed predictor, the no-GNRI forest is the same in every dataset: fit it once
  varying <- ms$pen_vars[vapply(ms$pen_vars, function(v) {
    any(vapply(imps[-1], function(d) !isTRUE(all.equal(d[[v]], imps[[1]][[v]])), logical(1)))
  }, logical(1))]
  single_ref <- identical(varying, "GNRI")
  no_gnri_vars <- setdiff(c(ms$pen_vars, ms$design_vars), "GNRI")
  no_gnri <- lapply(if (single_ref) ids[1] else ids, function(i) {
    p <- file.path(cache_dir, sprintf("noGNRI_cv_%d.rds", i))
    if (file.exists(p)) {
      return(readRDS(p))
    }
    r <- rsf_oof(imps[[i]], fold_id, ms, no_gnri_vars)
    saveRDS(r, p)
    r
  })
  no_gnri_for <- function(i) if (single_ref) no_gnri[[1]] else no_gnri[[i]]

  rows <- list()
  ref_cache <- list()
  for (i in ids) {
    for (cs in names(codes)) {
      for (g in names(groups)) {
        rows[[length(rows) + 1]] <- cbind(
          model = "full", imputation = i, cause = cs, group = g, months = cfg$horizons_months,
          score_subset(base, oof[[i]][[cs]], groups[[g]], codes[[cs]], hz)
        )
        key <- paste(if (single_ref) "one" else i, cs, g)
        if (is.null(ref_cache[[key]])) ref_cache[[key]] <- score_subset(base, no_gnri_for(i)[[cs]], groups[[g]], codes[[cs]], hz)
        rows[[length(rows) + 1]] <- cbind(
          model = "noGNRI", imputation = i, cause = cs, group = g, months = cfg$horizons_months, ref_cache[[key]]
        )
      }
    }
  }
  metrics <- do.call(rbind, rows)

  p2 <- read.csv(results_path("phase2_rsf_cv_by_dataset.csv"))
  chk <- merge(metrics[metrics$model == "full" & metrics$group == "all", c("imputation", "cause", "months", "AUC")],
    p2[, c("imputation", "cause", "months", "AUC")],
    by = c("imputation", "cause", "months"), suffixes = c("", "_p2")
  )
  if (max(abs(chk$AUC - chk$AUC_p2), na.rm = TRUE) > 1e-8) stop("Full-RSF scores do not reproduce Phase 2.")

  j <- merge(metrics[metrics$model == "full", ], metrics[metrics$model == "noGNRI", ],
    by = c("imputation", "cause", "group", "months"), suffixes = c("_full", "_noG")
  )
  j$dAUC <- j$AUC_full - j$AUC_noG
  j$dIPA <- j$IPA_full - j$IPA_noG
  keep <- c("imputation", "cause", "months", "dAUC", "dIPA")
  w_all <- j[j$group == "all", keep]
  dd <- merge(j[j$group == "imputed", keep], j[j$group == "observed", keep],
    by = c("imputation", "cause", "months"), suffixes = c("_imp", "_obs")
  )
  dd <- merge(dd, setNames(w_all, c("imputation", "cause", "months", "dAUC_all", "dIPA_all")),
    by = c("imputation", "cause", "months")
  )
  dd$DiD_AUC <- dd$dAUC_imp - dd$dAUC_obs
  dd$DiD_IPA <- dd$dIPA_imp - dd$dIPA_obs
  did <- do.call(rbind, lapply(split(dd, list(dd$cause, dd$months), drop = TRUE), function(z) {
    data.frame(
      cause = z$cause[1], months = z$months[1],
      dAUC_all = mean(z$dAUC_all), dAUC_observed = mean(z$dAUC_obs),
      dAUC_imputed = mean(z$dAUC_imp), DiD_AUC = mean(z$DiD_AUC),
      DiD_AUC_sd_across = sd(z$DiD_AUC),
      dIPA_observed = mean(z$dIPA_obs), dIPA_imputed = mean(z$dIPA_imp),
      DiD_IPA = mean(z$DiD_IPA)
    )
  }))
  did <- did[order(did$cause, did$months), ]

  # The across-dataset SD is imputation spread; this bootstrap gives the sampling uncertainty
  jh <- match(cfg$leak_boot_months, cfg$horizons_months)
  i1 <- ids[1]
  set.seed(cfg$seed + 1)
  auc1 <- function(risk, idx, cause) score_subset(base, risk[, jh, drop = FALSE], idx, cause, hz[jh])$AUC
  boot <- do.call(rbind, lapply(names(codes), function(cs) {
    rf <- oof[[i1]][[cs]]
    rn <- no_gnri_for(i1)[[cs]]
    did_b <- replicate(cfg$leak_boot_B, {
      b <- sample.int(n, n, replace = TRUE)
      bo <- b[!r_imp[b]]
      bi <- b[r_imp[b]]
      (auc1(rf, bi, codes[[cs]]) - auc1(rn, bi, codes[[cs]])) - (auc1(rf, bo, codes[[cs]]) - auc1(rn, bo, codes[[cs]]))
    })
    q <- quantile(did_b, c(0.025, 0.975), na.rm = TRUE)
    data.frame(
      cause = cs, months = cfg$leak_boot_months, dataset = i1,
      DiD_point = dd$DiD_AUC[dd$cause == cs & dd$months == cfg$leak_boot_months & dd$imputation == i1],
      boot_se = sd(did_b, na.rm = TRUE), lo95 = q[1], hi95 = q[2],
      n_boot_ok = sum(is.finite(did_b)), row.names = NULL
    )
  }))

  all_g <- unlist(lapply(imps, function(d) d$GNRI))
  covs <- setdiff(c(ms$pen_vars, ms$design_vars), "GNRI")
  cox <- do.call(rbind, lapply(names(codes), function(cs) {
    fm <- as.formula(paste0("Surv(cr_time, ", cause_codes()[[cs]]$event, ") ~ ", paste(covs, collapse = " + "), " + GNRI_s * Rimp"))
    e <- NULL
    v <- NULL
    for (i in seq_along(imps)) {
      d <- imps[[i]]
      d$GNRI_s <- (d$GNRI - mean(all_g)) / sd(all_g)
      d$Rimp <- as.integer(r_imp)
      f <- survival::coxph(fm, data = d, ties = "efron")
      b <- coef(f)
      vc <- vcov(f)
      t1 <- "GNRI_s"
      t2 <- "GNRI_s:Rimp"
      e <- rbind(e, c(
        GNRI_observed = unname(b[t1]), GNRI_imputed = unname(b[t1] + b[t2]),
        interaction = unname(b[t2]), imputed_indicator = unname(b["Rimp"])
      ))
      v <- rbind(v, c(vc[t1, t1], vc[t1, t1] + vc[t2, t2] + 2 * vc[t1, t2], vc[t2, t2], vc["Rimp", "Rimp"]))
    }
    colnames(v) <- colnames(e)
    pr <- rubin_pool(e, v, nu_com = n - length(coef(f)))
    pr <- data.frame(term = pr$variable, pr[, c("est", "se", "lo", "hi", "p", "fmi")])
    pr$HR <- exp(pr$est)
    pr$HR_lo <- exp(pr$lo)
    pr$HR_hi <- exp(pr$hi)
    cbind(cause = cs, pr)
  }))

  list(metrics = metrics, did = did, bootstrap = boot, cox = cox)
}
