# 00_run_all.R: run the full pipeline, each script in a fresh R session, and record the session.
# In: data/raw/   Out: data/derived/, results/, output/, output/session_info.txt

library(here)
source(here("R", "config.R"))

steps <- c("01_prepare_data.R", "02_fit_models.R", "03_diagnostics.R", "04_figures_tables.R")
for (s in steps) {
  message("Running ", s)
  if (system2("Rscript", here("scripts", s)) != 0) stop(s, " failed")
}

# Attach every package the pipeline uses so sessionInfo() lists their versions
pkgs <- c("survival", "mice", "rstan", "riskRegression", "prodlim", "randomForestSRC", "label.switching")
invisible(lapply(pkgs, library, character.only = TRUE))
out <- project_paths()$output
dir.create(out, showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
