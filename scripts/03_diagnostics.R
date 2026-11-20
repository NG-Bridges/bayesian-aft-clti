# 03_diagnostics.R: calibration, between-cause and mixture contrasts, mixture stability, imputation leakage, hazard form.
# In: data/derived/prepped.rds, results/ from 02   Out: results/{calibration,stacked_contrasts,mixture_ari,leakage_check,hazard_form}/

library(here)
suppressMessages(invisible(lapply(c("survival", "prodlim", "riskRegression", "randomForestSRC"), library, character.only = TRUE)))
for (f in c(
  "config.R", "model_setup.R", "checkpoint.R", "fit.R", "criteria.R", "predict.R", "rsf.R",
  "calibration.R", "contrasts.R", "mixture_stability.R", "leakage.R", "hazard_form.R"
)) {
  source(here("R", f))
}

cfg <- analysis_config()
set.seed(cfg$seed)
options(rf.cores = cfg$rsf_cores, mc.cores = 1)
dat <- load_prepped()
ms <- model_setup(dat, cfg)
fold_id <- make_folds(ms$base$cr_status, cfg$n_folds, cfg$seed)
ckpt_dir <- current_checkpoint_dir()
out <- function(...) {
  p <- results_path(...)
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  p
}
pooled <- read.csv(results_path("phase1_coefficients_pooled.csv"), stringsAsFactors = FALSE)

message("Calibration")
oof <- lapply(setNames(cfg$calib_models, cfg$calib_models), function(mdl) {
  lapply(seq_len(cfg$n_imputations), function(i) load_oof(mdl, i, ckpt_dir, cfg))
})
avail <- Reduce(intersect, lapply(oof, function(l) which(!vapply(l, is.null, logical(1)))))
if (length(avail) < cfg$n_imputations) stop("Out-of-fold predictions missing; run 02_fit_models.R --phases=2,3.")
write.csv(verify_oof(oof, avail, ms$imps, cfg), out("calibration", "verification.csv"), row.names = FALSE)
cal <- calibration_analysis(oof, pseudo_values(ms$base, cfg), avail, cfg)
write.csv(cal$by_dataset, out("calibration", "calibration_by_dataset.csv"), row.names = FALSE)
write.csv(cal$spread, out("calibration", "calibration_imputation_spread.csv"), row.names = FALSE)
write.csv(cal$summary, out("calibration", "calibration_summary.csv"), row.names = FALSE)
write.csv(cal$deciles, out("calibration", "calibration_deciles.csv"), row.names = FALSE)
write.csv(cal$curves, out("calibration", "calibration_curves.csv"), row.names = FALSE)
# Patient-level averaged predictions for the Figure 3 histograms; stays in results/
saveRDS(cal$pbar, out("calibration", "averaged_predictions.rds"))
rm(oof, cal)

message("Between-cause and mixture-class contrasts")
write.csv(between_cause_differences(ckpt_dir, pooled, cfg$n_imputations),
  out("stacked_contrasts", "between_cause_differences.csv"),
  row.names = FALSE
)
write.csv(mixture_class_differences(ckpt_dir, pooled, cfg$n_imputations),
  out("stacked_contrasts", "mixture_class_differences.csv"),
  row.names = FALSE
)

message("Mixture stability")
diag <- read.csv(results_path("phase1_fit_diagnostics.csv"))
ari <- lapply(setNames(names(cause_codes()), names(cause_codes())), mixture_stability, ckpt_dir = ckpt_dir, diag = diag, ms = ms)
for (cs in names(ari)) {
  write.csv(ari[[cs]]$pairs, out("mixture_ari", sprintf("ari_pairs_%s.csv", cs)), row.names = FALSE)
  write.csv(ari[[cs]]$stability, out("mixture_ari", sprintf("patient_stability_%s.csv", cs)), row.names = FALSE)
}
write.csv(do.call(rbind, lapply(ari, `[[`, "summary")), out("mixture_ari", "ari_summary.csv"), row.names = FALSE)

message("Imputation leakage")
lk <- leakage_check(ms, fold_id, ckpt_dir, dirname(out("leakage_check", "x")))
write.csv(lk$metrics, out("leakage_check", "leakage_metrics_by_dataset.csv"), row.names = FALSE)
write.csv(lk$did, out("leakage_check", "leakage_did_summary.csv"), row.names = FALSE)
write.csv(lk$bootstrap, out("leakage_check", "leakage_bootstrap.csv"), row.names = FALSE)
write.csv(lk$cox, out("leakage_check", "leakage_cox_interaction.csv"), row.names = FALSE)

message("Hazard form")
hf <- hazard_form_tests(dat$imps, dat$P$core, cfg)
write.csv(hf$pooled, out("hazard_form", "hazard_form_pooled.csv"), row.names = FALSE)
write.csv(hf$per_draw, out("hazard_form", "hazard_form_per_draw.csv"), row.names = FALSE)
