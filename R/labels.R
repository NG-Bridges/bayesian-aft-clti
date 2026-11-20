#' Display labels for the design indicators (0/1 as recoded in 02_fit_models.R).
design_labels <- function() {
  c(COHORT = "Cohort 2 (no adequate GSV)", ARM = "Endovascular arm", FACTOR_B = "Infrapopliteal disease")
}

#' Display labels for the clinical predictors, in table order.
#' The column ESKD is CKDGRADE 5, which identifies kidney transplant recipients (docs/decisions.md).
clinical_labels <- function() {
  c(
    Wound = "WIfI wound grade", ISCHEMIA = "WIfI ischemia grade",
    INFECTION = "WIfI foot-infection grade",
    PREV_INFR_IND = "Prior infrainguinal revascularization",
    AGE_T = "Age", GENDER = "Sex (female)", AMB_STAT3CAT = "Ambulatory status",
    SMOKE_HX = "Smoking history", DIAB_YN = "Diabetes",
    prior_cardiac_dx = "Prior cardiac disease", CHF_NYHA = "Heart failure (NYHA class)",
    cerebrovascular_dx = "Cerebrovascular disease", COPD_YN = "COPD",
    ESKD = "Kidney transplant", eGFR = "eGFR", eGFR_sq = "eGFR$^2$",
    GNRI = "GNRI", STATIN = "Statin use"
  )
}

model_labels <- function() c(AFT_single = "Single-component", AFT_frailty = "Frailty", AFT_mixture = "Gated mixture")

cause_labels <- function() c(amputation = "Amputation", death = "Death without amputation")
