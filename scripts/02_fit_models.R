# 02_fit_models.R: Bayesian AFT fits (WAIC/DIC3, pooled coefficients), RSF importance and CV, Bayesian CV.
# In: data/derived/prepped.rds   Out: results/checkpoints/, results/cv_predictions/, results/phase{1,2,3}_*, results/predictor_moments.csv

library(here)
suppressMessages(invisible(lapply(
  c("rstan", "survival", "riskRegression", "prodlim", "randomForestSRC", "label.switching"),
  library,
  character.only = TRUE
)))
for (f in c("config.R", "model_setup.R", "checkpoint.R", "fit.R", "criteria.R", "predict.R", "rsf.R")) {
  source(here("R", f))
}
rstan_options(auto_write = TRUE)

# Usage: Rscript scripts/02_fit_models.R [--phases=1,2,3] [--cv-imputations=50]
args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit) else default
}
run_phases <- as.integer(strsplit(arg_value("phases", "1,2,3"), ",")[[1]])

cfg <- analysis_config()
cfg$cv_n_imputations <- as.integer(arg_value("cv-imputations", cfg$cv_n_imputations))
set.seed(cfg$seed)
options(rf.cores = cfg$rsf_cores, mc.cores = 1) # no forking: the original runs were on Windows

out_dir <- project_paths()$results
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
causes <- cause_codes()
dat <- load_prepped()
if (max(cfg$n_imputations, cfg$cv_n_imputations, cfg$rsf_n_imputations) > length(dat$imps)) {
  stop("Requested more imputations than 01_prepare_data.R produced (", length(dat$imps), ").")
}
ms <- model_setup(dat, cfg)
base <- ms$base
imps <- ms$imps
rm(dat)

write.csv(data.frame(variable = ms$pen_vars, mean = as.numeric(ms$fix_mom$mu), sd = as.numeric(ms$fix_mom$sd)),
  file.path(out_dir, "predictor_moments.csv"),
  row.names = FALSE
)
fold_id <- make_folds(base$cr_status, cfg$n_folds, cfg$seed)
ckpt_dir <- checkpoint_dir(ms)
ckpt <- function(...) checkpoint_path(ckpt_dir, ...)
message("Checkpoints: ", ckpt_dir)

# A systematic scoring failure would otherwise surface as NA columns after days of fitting
local({
  set.seed(cfg$seed)
  rr <- matrix(runif(nrow(base) * length(cfg$horizons)), nrow(base))
  st <- score_cr(base$cr_time, base$cr_status, causes$amputation$code, rr, cfg$horizons)
  if (all(is.na(st$AUC))) stop("Scoring self-test failed: riskRegression::Score() returned nothing.")
})

inference_data <- function(i) prepare_design(imps[[i]], ms, ms$fix_mom, cfg)

if (1 %in% run_phases) {
  message("Phase 1: ", length(cfg$models) * length(causes) * cfg$n_imputations, " full-sample fits")
  ic_rows <- list()
  diag_rows <- list()
  elpd_store <- list()
  draw_store <- list()
  for (mdl in cfg$models) {
    for (cs in names(causes)) {
      ecol <- causes[[cs]]$event
      for (i in seq_len(cfg$n_imputations)) {
        r <- checkpointed(ckpt("inf", mdl, cs, i), function() {
          d <- inference_data(i)
          fit <- run_stan(mdl, d, ecol, cs, ms)
          summarise_inference_fit(fit, mdl, d, ecol, ms)
        })
        ic_rows[[length(ic_rows) + 1]] <- data.frame(
          model = mdl, cause = cs, imputation = i,
          waic = r$ic$waic, se_waic = r$ic$se_waic, p_waic = r$ic$p_waic,
          dic3 = r$ic$dic3, p_dic3 = r$ic$p_dic3, lppd = r$ic$lppd, Dbar = r$ic$Dbar,
          max_p_waic_i = r$ic$max_p_waic_i, n_p_waic_gt_04 = r$ic$n_p_waic_gt_04
        )
        diag_rows[[length(diag_rows) + 1]] <- cbind(model = mdl, cause = cs, imputation = i, r$diag)
        elpd_store[[paste(mdl, cs, i)]] <- r$ic$elpd_i
        draw_store[[paste(mdl, cs, i)]] <- r$draws
      }
    }
  }
  ic_tab <- do.call(rbind, ic_rows)
  diag_tab <- do.call(rbind, diag_rows)
  write.csv(ic_tab, file.path(out_dir, "phase1_information_criteria.csv"), row.names = FALSE)
  write.csv(diag_tab, file.path(out_dir, "phase1_fit_diagnostics.csv"), row.names = FALSE)
  if (any(!diag_tab$converged)) {
    message(sum(!diag_tab$converged), "/", nrow(diag_tab), " fits fail the convergence gate; see phase1_fit_diagnostics.csv")
  }

  # Frailty quadrature error at the fitted posterior, against 4x the nodes (2x is not yet
  # converged at the amputation sigma_u). Single draws can be off by several deviance units;
  # WAIC is what the comparison uses, and it moves far less.
  d1 <- inference_data(1)
  xp1 <- design_mat(d1, ms$model_pen_vars)
  xd1 <- design_mat(d1, ms$design_vars)
  for (cs in names(causes)) {
    m <- draw_store[[paste("AFT_frailty", cs, 1)]]
    ev1 <- as.integer(d1[[causes[[cs]]$event]])
    lp <- m[, "alpha"] + m[, ms$model_pen_vars] %*% t(xp1) + m[, ms$design_vars] %*% t(xd1)
    ll <- lapply(c(1, 4) * cfg$gh_nodes, function(nodes) {
      t(vapply(seq_len(nrow(m)), function(s) {
        gh_loglik(log(d1$cr_time), ev1, lp[s, ], m[s, "k_shape"], m[s, "sigma_u"], nodes)
      }, numeric(nrow(d1))))
    })
    message(sprintf(
      "Quadrature check (%s, %d vs %d nodes): WAIC shift %.2f, largest per-draw deviance shift %.2f",
      cs, cfg$gh_nodes, 4 * cfg$gh_nodes,
      abs(ic_from_loglik(ll[[1]])$waic - ic_from_loglik(ll[[2]])$waic),
      2 * max(abs(rowSums(ll[[1]]) - rowSums(ll[[2]])))
    ))
  }
  rm(d1, xp1, xd1, lp, ll)

  icsum <- do.call(rbind, lapply(split(ic_tab, list(ic_tab$model, ic_tab$cause), drop = TRUE), function(z) {
    data.frame(
      model = z$model[1], cause = z$cause[1],
      waic = mean(z$waic), waic_sd = sd(z$waic),
      dic3 = mean(z$dic3), dic3_sd = sd(z$dic3),
      p_waic = mean(z$p_waic), p_dic3 = mean(z$p_dic3),
      mean_se_waic = mean(z$se_waic),
      max_p_waic_i = max(z$max_p_waic_i), n_flagged_obs = sum(z$n_p_waic_gt_04)
    )
  }))
  icsum <- icsum[order(icsum$cause, icsum$waic), ]
  write.csv(icsum, file.path(out_dir, "phase1_ic_summary.csv"), row.names = FALSE)

  # Paired within-dataset differences: every model sees identical data within a dataset
  pair_rows <- list()
  for (cs in names(causes)) {
    for (cb in combn(cfg$models, 2, simplify = FALSE)) {
      a <- cb[1]
      b <- cb[2]
      dw <- numeric(0)
      dd <- numeric(0)
      within_se <- numeric(0)
      for (i in seq_len(cfg$n_imputations)) {
        za <- ic_tab[ic_tab$model == a & ic_tab$cause == cs & ic_tab$imputation == i, ]
        zb <- ic_tab[ic_tab$model == b & ic_tab$cause == cs & ic_tab$imputation == i, ]
        dw <- c(dw, za$waic - zb$waic)
        dd <- c(dd, za$dic3 - zb$dic3)
        ea <- elpd_store[[paste(a, cs, i)]]
        eb <- elpd_store[[paste(b, cs, i)]]
        within_se <- c(within_se, 2 * sqrt(length(ea)) * sd(ea - eb))
      }
      pair_rows[[length(pair_rows) + 1]] <- data.frame(
        cause = cs, comparison = paste(a, "-", b),
        d_waic = mean(dw), d_waic_sd_across = sd(dw), d_waic_se_within = mean(within_se),
        d_dic3 = mean(dd), d_dic3_sd_across = sd(dd),
        favours = if (mean(dw) < -2 * mean(within_se)) a else if (mean(dw) > 2 * mean(within_se)) b else "tie"
      )
    }
  }
  pair_tab <- do.call(rbind, pair_rows)
  write.csv(pair_tab, file.path(out_dir, "phase1_ic_paired_differences.csv"), row.names = FALSE)

  # Stacking draws across imputations is exact for the MI posterior; Rubin's rules only supply FMI
  coef_rows <- list()
  rubin_rows <- list()
  for (mdl in cfg$models) {
    for (cs in names(causes)) {
      keys <- paste(mdl, cs, seq_len(cfg$n_imputations))
      m_all <- do.call(rbind, draw_store[keys])
      q <- t(apply(m_all, 2, quantile, c(0.5, 0.025, 0.975), names = FALSE))
      stem <- sub("_c[0-9]+$", "", colnames(m_all))
      coef_rows[[length(coef_rows) + 1]] <- data.frame(
        model = mdl, cause = cs, variable = colnames(m_all),
        block = ifelse(stem %in% ms$design_vars, "unpenalized", ifelse(stem %in% ms$model_pen_vars, "penalized", "ancillary")),
        post_mean = colMeans(m_all), post_median = q[, 1], lo95 = q[, 2], hi95 = q[, 3],
        p_direction = apply(m_all, 2, function(x) max(mean(x > 0), mean(x < 0))), row.names = NULL
      )
      e <- t(sapply(keys, function(k) colMeans(draw_store[[k]])))
      v <- t(sapply(keys, function(k) apply(draw_store[[k]], 2, var)))
      rubin_rows[[length(rubin_rows) + 1]] <- cbind(model = mdl, cause = cs, rubin_pool(e, v, nrow(base) - ncol(e)))
    }
  }
  coef_tab <- do.call(rbind, coef_rows)
  rubin_tab <- do.call(rbind, rubin_rows)
  write.csv(coef_tab, file.path(out_dir, "phase1_coefficients_pooled.csv"), row.names = FALSE)
  write.csv(rubin_tab, file.path(out_dir, "phase1_coefficients_rubin.csv"), row.names = FALSE)
  saveRDS(
    list(ic = ic_tab, diag = diag_tab, coef = coef_tab, rubin = rubin_tab, pairs = pair_tab),
    file.path(out_dir, "phase1_results.rds")
  )
  rm(elpd_store, draw_store)
}

if (2 %in% run_phases) {
  message("Phase 2: RSF importance and 5-fold CV on ", cfg$rsf_n_imputations, " datasets")
  imp_rows <- list()
  for (i in seq_len(cfg$rsf_n_imputations)) {
    vi <- checkpointed(ckpt("rsf_imp", i), function() rsf_importance(imps[[i]], ms))
    imp_rows[[length(imp_rows) + 1]] <- cbind(imputation = i, vi$tab)
  }
  rsf_imp <- do.call(rbind, imp_rows)
  rsf_imp_summ <- do.call(rbind, lapply(split(rsf_imp, rsf_imp$variable), function(z) {
    data.frame(
      variable = z$variable[1],
      amp_vimp = mean(z$amputation, na.rm = TRUE), amp_sd = sd(z$amputation, na.rm = TRUE),
      death_vimp = mean(z$death, na.rm = TRUE), death_sd = sd(z$death, na.rm = TRUE)
    )
  }))
  rsf_imp_summ$amp_rank <- rank(-rsf_imp_summ$amp_vimp)
  rsf_imp_summ$death_rank <- rank(-rsf_imp_summ$death_vimp)
  write.csv(rsf_imp, file.path(out_dir, "phase2_rsf_importance_by_dataset.csv"), row.names = FALSE)
  write.csv(rsf_imp_summ, file.path(out_dir, "phase2_rsf_importance_pooled.csv"), row.names = FALSE)

  rsf_rows <- list()
  for (i in seq_len(cfg$rsf_n_imputations)) {
    oof <- checkpointed(ckpt("rsf_cv", i), function() rsf_oof(imps[[i]], fold_id, ms))
    d <- imps[[i]]
    for (cs in names(causes)) {
      sc <- score_cr(d$cr_time, d$cr_status, causes[[cs]]$code, oof[[cs]], cfg$horizons)
      rsf_rows[[length(rsf_rows) + 1]] <- cbind(model = "RSF", cause = cs, imputation = i, months = cfg$horizons_months, sc)
    }
  }
  rsf_cv <- do.call(rbind, rsf_rows)
  write.csv(rsf_cv, file.path(out_dir, "phase2_rsf_cv_by_dataset.csv"), row.names = FALSE)
  saveRDS(rsf_cv, file.path(out_dir, "phase2_rsf_cv.rds"))
}

if (3 %in% run_phases) {
  message(
    "Phase 3: ", length(cfg$cv_models) * length(causes) * cfg$cv_n_imputations * cfg$n_folds,
    " cross-validation fits; out-of-fold CIFs are scored once per dataset"
  )
  # Out-of-fold CIFs for one dataset, or NULL if a fold's fit failed
  cv_oof <- function(mdl, i) {
    d <- imps[[i]]
    a <- matrix(NA_real_, nrow(d), length(cfg$horizons))
    dm <- a
    for (k in seq_len(cfg$n_folds)) {
      tr <- d[fold_id != k, ]
      te <- d[fold_id == k, ]
      mom <- std_moments(tr, ms$pen_vars)
      tr_s <- prepare_design(tr, ms, mom, cfg)
      te_s <- prepare_design(te, ms, mom, cfg)
      drs <- list()
      for (cs in names(causes)) {
        # One failed fit skips its dataset rather than ending a multi-day run
        drs[[cs]] <- tryCatch(
          checkpointed(ckpt("cv", mdl, cs, i, "f", k), function() {
            fit <- run_stan(mdl, tr_s, causes[[cs]]$event, cs, ms)
            out <- prediction_draws(fit, mdl, cfg$pred_draws)
            out$diag <- fit_diagnostics(fit, mdl, cfg)
            out
          }),
          missing_checkpoint = function(e) stop(e),
          error = function(e) {
            message(sprintf("Fit failed: %s/%s/imputation %d/fold %d: %s", mdl, cs, i, k, conditionMessage(e)))
            NULL
          }
        )
        if (is.null(drs[[cs]])) {
          return(NULL)
        }
      }
      pr <- predict_cif(
        drs$amputation, drs$death,
        design_mat(te_s, ms$model_pen_vars), design_mat(te_s, ms$design_vars), design_mat(te_s, ms$gate_vars), ms
      )
      a[fold_id == k, ] <- pr$amputation
      dm[fold_id == k, ] <- pr$death
    }
    list(amputation = a, death = dm)
  }

  # Saved per dataset (patient-level, gitignored) so a rerun resumes and 03_diagnostics.R reuses them
  pred_dir <- file.path(out_dir, "cv_predictions", basename(ckpt_dir))
  dir.create(pred_dir, showWarnings = FALSE, recursive = TRUE)
  cv_rows <- list()
  for (mdl in cfg$cv_models) {
    for (i in seq_len(cfg$cv_n_imputations)) {
      f <- file.path(pred_dir, sprintf("oof_%s_grid%d_%d.rds", mdl, cfg$grid_n, i))
      if (file.exists(f)) {
        oof <- readRDS(f)
      } else {
        oof <- cv_oof(mdl, i)
        if (is.null(oof)) {
          message("Skipping ", mdl, " dataset ", i, " (a fold failed)")
          next
        }
        saveRDS(oof, f)
      }
      d <- imps[[i]]
      for (cs in names(causes)) {
        sc <- score_cr(d$cr_time, d$cr_status, causes[[cs]]$code, oof[[cs]], cfg$horizons)
        cv_rows[[length(cv_rows) + 1]] <- cbind(model = mdl, cause = cs, imputation = i, months = cfg$horizons_months, sc)
      }
      message(sprintf("Phase 3: %s dataset %d/%d scored (%s)", mdl, i, cfg$cv_n_imputations, format(Sys.time(), "%H:%M")))
    }
  }
  bayes_cv <- do.call(rbind, cv_rows)
  write.csv(bayes_cv, file.path(out_dir, "phase3_bayes_cv_by_dataset.csv"), row.names = FALSE)
  saveRDS(bayes_cv, file.path(out_dir, "phase3_bayes_cv.rds"))
}
