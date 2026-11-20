# 01_prepare_data.R: assemble and clean the trial export, build the cohort, impute albumin, weight and height (m = 50).
# In: data/raw/sas/*.sas7bdat (example: simulated)   Out: data/derived/ export, CLI_clean.csv, manifest, prepped.rds

library(here)
for (f in c("config.R", "export.R", "example_data.R", "clean.R", "data.R")) source(here("R", f))

cfg <- analysis_config()
paths <- project_paths()
seed <- cfg$seed
m <- cfg$n_imputations
eps <- 0.5 / 365.25 # half a day, in years
# The published imputations were drawn in a session using L'Ecuyer-CMRG (recorded in the saved mids object)
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)

derived <- paths$derived
dir.create(derived, showWarnings = FALSE, recursive = TRUE)
if (example_mode()) {
  # A simulated export with the trial's layout and coding stands in for the SAS datasets
  sim <- simulate_trial_export(n = 400, seed = seed)
  write.csv(sim$export, file.path(derived, "FullAugmentedData.csv"), row.names = FALSE)
  write.csv(sim$randomization, file.path(derived, "Randomization.csv"), row.names = FALSE)
  rm(sim)
} else {
  export <- build_trial_export(file.path(paths$raw, "sas"))
  readr::write_csv(export$data, file.path(derived, "FullAugmentedData.csv"), na = "NA")
  readr::write_csv(export$randomization, file.path(derived, "Randomization.csv"))
  rm(export)
}
# Read back rather than reuse the data frames: the cleaning was built on read.csv's column typing
raw <- read.csv(file.path(derived, "FullAugmentedData.csv"), check.names = FALSE, stringsAsFactors = FALSE)
randomization <- read.csv(file.path(derived, "Randomization.csv"), check.names = FALSE, stringsAsFactors = FALSE)
message(sprintf("Trial export: %d records x %d columns", nrow(raw), ncol(raw)))

cleaned <- clean_trial_data(raw, randomization)
clean_path <- file.path(derived, "CLI_clean.csv")
write.csv(cleaned$data, clean_path, row.names = FALSE)
write.csv(cleaned$manifest, file.path(derived, "cleaning_manifest.csv"), row.names = FALSE)
message(sprintf("Cleaned data: %d patients x %d columns", nrow(cleaned$data), ncol(cleaned$data)))

# Read back rather than reuse cleaned$data: read.csv's integer/double typing is part of the data the models saw
ds <- analysis_cohort(cli_load(clean_path))
message(sprintf("Analysis cohort: n = %d", nrow(ds)))

# 12 patients have cr_time == 0, where a log-time model is undefined
for (v in c("cr_time", "amp_time", "death_time")) ds[[v]] <- pmax(ds[[v]], eps)

ds$ev_amp <- as.integer(ds$cr_status == 1)
ds$ev_death <- as.integer(ds$cr_status == 2)
ds$ev_reint <- ds$EVNT_REINTIOS

# Event indicators and cumulative hazards keep the imputation model congenial
# with the survival analysis (White and Royston 2009)
ds$NA_amp <- nelson_aalen(ds$cr_time, ds$ev_amp)
ds$NA_death <- nelson_aalen(ds$cr_time, ds$ev_death)
ds$NA_reint <- nelson_aalen(ds$TMTO_REINTIOS, ds$ev_reint)

targets <- c("ALBUMIN", "WEIGHT", "HEIGHT")
auxiliary <- c("ev_amp", "ev_death", "ev_reint", "NA_amp", "NA_death", "NA_reint")
mi_in <- ds[, c(targets, core_predictors(), auxiliary)]
message(sprintf(
  "Missing: albumin %d (%.1f%%), weight %d, height %d",
  sum(is.na(ds$ALBUMIN)), 100 * mean(is.na(ds$ALBUMIN)), sum(is.na(ds$WEIGHT)), sum(is.na(ds$HEIGHT))
))

predictors <- matrix(0, ncol(mi_in), ncol(mi_in), dimnames = list(names(mi_in), names(mi_in)))
predictors[targets, ] <- 1
diag(predictors) <- 0
method <- setNames(rep("", ncol(mi_in)), names(mi_in))
method[targets] <- "pmm"

imp <- mice::mice(mi_in,
  m = m, maxit = 5, method = method, predictorMatrix = predictors,
  seed = seed, printFlag = FALSE
)

imputed <- lapply(seq_len(m), function(i) {
  di <- ds
  ci <- mice::complete(imp, i)
  di[targets] <- ci[targets]
  di$GNRI <- gnri(di$ALBUMIN, di$WEIGHT, di$HEIGHT, di$GENDER)
  di
})
if (anyNA(unlist(lapply(imputed, `[[`, "GNRI")))) stop("GNRI is missing after imputation")

# base_ds keeps the unimputed values, so downstream scripts can tell measured from imputed patients
saveRDS(
  list(P = list(core = core_predictors()), base_ds = ds, imputed = imputed, imp_obj = imp),
  file.path(derived, "prepped.rds")
)
message("Wrote ", file.path(derived, "prepped.rds"))
