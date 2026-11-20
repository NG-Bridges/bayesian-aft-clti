# 04_figures_tables.R: the article's figures and the LaTeX tables, from saved summaries and Phase 1 draws.
# In: results/ from 02 and 03 (summaries, calibration, contrasts, Phase 1 draws)   Out: output/figures/, output/tables/

library(here)
for (f in c("config.R", "labels.R", "latex.R", "tables.R", "figures.R", "fit.R", "checkpoint.R")) source(here("R", f))

cfg <- analysis_config()
read_result <- function(...) {
  p <- results_path(...)
  if (!file.exists(p)) stop("Missing ", p)
  read.csv(p, stringsAsFactors = FALSE)
}
fig_dir <- file.path(project_paths()$output, "figures")
tab_dir <- file.path(project_paths()$output, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

co <- read_result("phase1_coefficients_pooled.csv")
draws <- stacked_single_draws(current_checkpoint_dir(), co[co$model == "AFT_single", ], cfg$n_imputations)
figure_aft_forest(
  co[co$model == "AFT_single", ], read_result("predictor_moments.csv"),
  read_result("phase2_rsf_importance_pooled.csv"), draws, file.path(fig_dir, "figure_aft_forest")
)
rm(draws)

# File name matches \includegraphics in the article
figure_cv(
  read_result("phase3_bayes_cv_by_dataset.csv"), read_result("phase2_rsf_cv_by_dataset.csv"),
  file.path(fig_dir, "figure_cv_AFT_vs_rsf")
)

for (mdl in names(model_labels())) {
  for (cs in names(cause_labels())) {
    v <- c(names(design_labels()), names(clinical_labels()))
    if (mdl == "AFT_mixture") v <- c(paste0(v, "_c1"), paste0(v, "_c2"))
    invisible(lapply(v, function(x) coef_row(co, mdl, cs, x)))
  }
}

ic <- ic_values(read_result("phase1_ic_summary.csv"), read_result("phase1_ic_paired_differences.csv"), cfg$n_imputations)
write.csv(ic, results_path("table_ic_values.csv"), row.names = FALSE)
writeLines(table_ic_tex(ic, cfg$n_imputations), file.path(tab_dir, "table_ic.tex"))

writeLines(
  table_coef_mixture_tex(co, read_result("stacked_contrasts", "mixture_class_differences.csv")),
  file.path(tab_dir, "tableS_coef_mixture.tex")
)
writeLines(table_coef_single_frailty_tex(co), file.path(tab_dir, "tableS_coef_single_frailty.tex"))

cv <- convergence_values(read_result("phase1_fit_diagnostics.csv"), cfg)
write.csv(
  do.call(rbind, lapply(names(cv), function(cs) data.frame(cause = cs, metric = rownames(cv[[cs]]), cv[[cs]]))),
  results_path("table_convergence_values.csv"),
  row.names = FALSE
)
writeLines(table_convergence_tex(cv, cfg), file.path(tab_dir, "table_convergence.tex"))

writeLines(
  table_swings_tex(read_result("stacked_contrasts", "between_cause_differences.csv")),
  file.path(tab_dir, "table_swings.tex")
)

cal <- lapply(c(summ = "calibration_summary.csv", dec = "calibration_deciles.csv", crv = "calibration_curves.csv"), function(f) read_result("calibration", f))
pbar_file <- results_path("calibration", "averaged_predictions.rds")
if (!file.exists(pbar_file)) stop("Missing ", pbar_file, "; run scripts/03_diagnostics.R first.")
pbar <- readRDS(pbar_file)
figure_calibration(cal$summ, cal$dec, cal$crv, pbar, cfg$calib_main_months,
  file.path(fig_dir, "figure_calibration_main.png"),
  width = 560 * length(cfg$calib_main_months) + 80, height = 1180, cfg = cfg
)
figure_calibration(cal$summ, cal$dec, cal$crv, pbar, cfg$horizons_months,
  file.path(fig_dir, "figureS_calibration_all.png"),
  width = 520 * length(cfg$horizons_months) + 80, height = 1100, cfg = cfg
)
