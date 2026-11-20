#' Clinical and design predictors shared by every model; GNRI and COHORT are added downstream.
core_predictors <- function() {
  # ESKD is CKDGRADE 5, which identifies kidney transplant recipients (docs/decisions.md).
  # The name is kept because the saved model fits are keyed on predictor names.
  c(
    "Wound", "ISCHEMIA", "INFECTION", "AGE_T", "GENDER", "eGFR", "ESKD",
    "prior_cardiac_dx", "CHF_NYHA", "DIAB_YN", "AMB_STAT3CAT", "COPD_YN",
    "cerebrovascular_dx", "SMOKE_HX", "STATIN", "PREV_INFR_IND", "ARM", "FACTOR_B"
  )
}

#' Geriatric Nutritional Risk Index (Bouillanne et al. 2005).
#' @param albumin Serum albumin, g/dL.
#' @param gender 1 = male, 2 = female.
gnri <- function(albumin, weight, height, gender) {
  # Lorentz ideal body weight, as in the original GNRI definition
  ideal <- ifelse(gender == 1,
    height - 100 - (height - 150) / 4,
    height - 100 - (height - 150) / 2.5
  )
  1.489 * (albumin * 10) + 41.7 * pmin(weight / ideal, 1)
}

#' Read the cleaned trial data and derive the composite predictors.
cli_load <- function(path) {
  d <- read.csv(path, check.names = FALSE)
  any_of <- function(cols) as.integer(do.call(pmax, c(d[, cols], list(na.rm = TRUE))) > 0)
  d$prior_cardiac_dx <- any_of(c("CAD_YN", "ANGINA_YN", "PREV_MI", "PREV_PCI", "PREV_CABG"))
  d$cerebrovascular_dx <- any_of(c("TIA_YN", "STROKE_YN"))
  d$GNRI <- gnri(d$ALBUMIN, d$WEIGHT, d$HEIGHT, d$GENDER)
  d$AGE_T <- (d$AGE - 40) / 10 # decades above 40
  d
}

#' Patients with every core predictor and a complete competing-risks outcome.
analysis_cohort <- function(d) {
  d[complete.cases(d[, c(core_predictors(), "cr_time", "cr_status")]), , drop = FALSE]
}

#' Marginal Nelson-Aalen cumulative hazard, evaluated at each subject's own time.
nelson_aalen <- function(time, event) {
  fit <- survival::survfit(survival::Surv(time, event) ~ 1)
  approx(fit$time, -log(pmax(fit$surv, 1e-12)),
    xout = time, method = "constant", rule = 2, f = 0
  )$y
}
