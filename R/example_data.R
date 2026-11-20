#' Parameters of the simulated example export.
#' Marked values are from the published BEST-CLI Table 1 and results (Farber et al., NEJM 2022);
#' the rest are round, plausible values. None is estimated from the trial data.
example_parameters <- function() {
  list(
    p_cohort2 = 396 / 1830, # Table 1
    age_mean = c(66.9, 68.6), age_sd = c(9.9, 9.2), # Table 1, by cohort
    p_female = 0.285, # Table 1
    bmi_mean = c(28.2, 27.0), bmi_sd = c(6, 5.7), # Table 1, by cohort
    p_diabetes = c(0.72, 0.60), # Table 1, by cohort
    p_current_smoker = 0.36, # Table 1
    p_dialysis = 0.11, # Table 1 (end-stage kidney disease)
    p_transplant = 0.012,
    p_statin = 0.70, p_aspirin = 0.67, p_clopidogrel = 0.22, # Table 1
    abi_mean = 0.58, # Table 1
    followup_mean = c(2.7, 1.6), # years, by cohort
    p_infrapopliteal = 0.66,
    albumin_mean = 3.6, albumin_sd = 0.5, p_albumin_missing = 0.46,
    # True AFT coefficients on the log-time scale (positive lengthens time to the event)
    beta_amp = c(
      intercept = 4.8, wound = -0.35, ischemia = -0.25, infection = -0.2, diabetes = -0.3,
      egfr10 = 0.05, gnri_sd = 0.25, endovascular = -0.3, infrapopliteal = -0.3, age10 = 0.1
    ),
    beta_death = c(
      intercept = 2.0, age10 = -0.35, chf = -0.4, copd = -0.3, ambulation = -0.3,
      egfr10 = 0.08, dialysis = -0.8, gnri_sd = 0.3, wound = -0.1
    ),
    shape = c(amputation = 0.7, death = 1.1)
  )
}

#' Simulate the two raw files, with the layout and coding of the BEST-CLI export.
#' Includes one case of each irregularity the cleaning handles (units, .x/.y pairs, near-empty forms,
#' patients with two procedure records).
#' @return List: `export` (as FullAugmentedData.csv) and `randomization` (as Randomization.csv).
simulate_trial_export <- function(n = 400, seed = 2026) {
  withr::with_seed(seed,
    .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
    .rng_sample_kind = "Rejection", simulate_export_draws(n, example_parameters())
  )
}

simulate_export_draws <- function(n, p) {
  bin <- function(prob, na = 0.01) with_na(rbinom(n, 1, prob), na)
  codes <- function(values, prob = NULL, na = 0.01) with_na(sample(values, n, TRUE, prob), na)
  num <- function(lo, hi, na = 0.03, digits = 1) with_na(round(runif(n, lo, hi), digits), na)
  with_na <- function(x, na) {
    x[runif(n) < na] <- NA
    x
  }
  clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

  cohort <- 1 + rbinom(n, 1, p$p_cohort2)
  arm <- sample(rep(1:2, length.out = n)) # raw 1 = surgical bypass, 2 = endovascular
  open <- arm == 1
  female <- rbinom(n, 1, p$p_female) == 1
  age <- round(clip(rnorm(n, p$age_mean[cohort], p$age_sd[cohort]), 40, 95), 1)
  height <- clip(rnorm(n, ifelse(female, 162, 176), 7), 145, 200)
  bmi <- clip(rnorm(n, p$bmi_mean[cohort], p$bmi_sd[cohort]), 16, 50)
  weight <- bmi * (height / 100)^2
  diabetes <- rbinom(n, 1, p$p_diabetes[cohort])
  renal <- sample(
    c("none", "dialysis", "transplant"), n, TRUE,
    c(1 - p$p_dialysis - p$p_transplant, p$p_dialysis, p$p_transplant)
  )
  egfr <- ifelse(renal == "dialysis", rnorm(n, 9, 3), ifelse(renal == "transplant", rnorm(n, 66, 15), rnorm(n, 70, 25)))
  egfr <- clip(egfr, 3, 150)
  ckd <- ifelse(renal == "dialysis", 4, ifelse(renal == "transplant", 5,
    cut(egfr, c(-Inf, 30, 45, 60, Inf), labels = FALSE) - 1
  ))
  wound <- sample(0:3, n, TRUE, c(0.25, 0.35, 0.25, 0.15))
  ischemia <- sample(1:4, n, TRUE, c(0.1, 0.2, 0.3, 0.4)) # exported 1-4
  infection <- sample(0:3, n, TRUE, c(0.6, 0.2, 0.15, 0.05))
  chf <- rbinom(n, 1, 0.15)
  copd <- rbinom(n, 1, 0.2)
  ambulation <- sample(1:3, n, TRUE, c(0.6, 0.3, 0.1))
  infrapopliteal <- rbinom(n, 1, p$p_infrapopliteal)
  albumin <- round(clip(rnorm(n, p$albumin_mean, p$albumin_sd), 1.5, 5), 1)
  gnri_true <- gnri(albumin, weight, height, ifelse(female, 2, 1))

  z_gnri <- (gnri_true - mean(gnri_true)) / sd(gnri_true)
  b <- p$beta_amp
  lp_amp <- b[["intercept"]] + b[["wound"]] * wound + b[["ischemia"]] * (ischemia - 1) + b[["infection"]] * infection +
    b[["diabetes"]] * diabetes + b[["egfr10"]] * (pmin(pmax(egfr, 5), 120) - 60) / 10 + b[["gnri_sd"]] * z_gnri +
    b[["endovascular"]] * (arm == 2) + b[["infrapopliteal"]] * infrapopliteal + b[["age10"]] * (age - 67) / 10
  b <- p$beta_death
  lp_death <- b[["intercept"]] + b[["age10"]] * (age - 67) / 10 + b[["chf"]] * chf + b[["copd"]] * copd +
    b[["ambulation"]] * (ambulation - 1) + b[["egfr10"]] * (pmin(pmax(egfr, 5), 120) - 60) / 10 +
    b[["dialysis"]] * (renal == "dialysis") + b[["gnri_sd"]] * z_gnri + b[["wound"]] * wound
  # Weibull AFT: log T = lp + W / shape, with W standard minimum extreme value
  t_amp <- exp(lp_amp) * rexp(n)^(1 / p$shape[["amputation"]])
  t_death <- exp(lp_death) * rexp(n)^(1 / p$shape[["death"]])
  censor <- runif(n, 0.1, 2 * p$followup_mean[cohort] - 0.1)
  t_reint <- rexp(n, ifelse(open, 0.08, 0.2))

  amp_event <- as.integer(t_amp < pmin(t_death, censor))
  amp_time <- ifelse(amp_event == 1, t_amp, pmin(t_death, censor))
  death_event <- as.integer(t_death < censor)
  death_time <- pmin(t_death, censor)
  reint_event <- as.integer(t_reint < pmin(t_death, censor))
  reint_time <- pmin(t_reint, t_death, censor)
  male_time <- pmin(ifelse(amp_event == 1, t_amp, Inf), t_reint, t_death, censor)
  male_event <- as.integer(male_time < censor)
  # Same-day events and censoring, which the analysis floors at half a day
  zero <- sample(n, 3)
  amp_time[zero] <- 0
  death_time[zero] <- 0

  height_in <- runif(n) < 0.3
  weight_lb <- runif(n) < 0.4
  creat <- round(clip(exp(log(1.1) - 0.9 * log(egfr / 70)), 0.3, 12), 2)
  creat_umol <- runif(n) < 0.1
  hba1c <- round(clip(ifelse(diabetes == 1, rnorm(n, 7.5, 1.5), rnorm(n, 5.7, 0.4)), 4, 14), 1)
  ifcc <- runif(n) < 0.1
  albumin_obs <- with_na(albumin, p$p_albumin_missing)
  both_sides <- runif(n) < 0.1
  wifi <- with_na(wound, 0.04)

  d <- data.frame(Random_Subnum = 900000 + seq_len(n))
  d$ARM <- arm
  d$ASSIGNED <- with_na(ifelse(runif(n) < 0.02, 3, arm), 0.01)
  d$RAND_ASSIGN <- arm
  d$ACTTRT <- with_na(arm, 0.03)
  d$STRATUM <- 2 * (cohort - 1) + 1 + infrapopliteal
  d$COHORT <- cohort
  d$SSGSV_YN <- as.integer(cohort == 1)
  d$AGE <- age
  d$GENDER <- ifelse(female, 2, 1)
  d$HISPANIC <- rbinom(n, 1, 0.1)
  d$BMI <- with_na(round(bmi, 1), 0.05)
  d$HEIGHT <- with_na(round(ifelse(height_in, height / 2.54, height), 1), 0.04)
  d$HEIGHT[sample(n, 2)] <- 0
  d$HEIGHT_UN <- ifelse(height_in, 2, 1)
  d$WEIGHT <- with_na(round(ifelse(weight_lb, weight / 0.45359237, weight), 1), 0.03)
  d$WEIGHT_UN <- ifelse(weight_lb, 2, 1)
  d$AMB_STAT <- codes(1:4)
  d$AMB_STAT3CAT <- with_na(ambulation, 0.01)
  d$LIVE_STAT <- codes(1:5, c(0.6, 0.2, 0.1, 0.05, 0.05))
  d$vascu_qol_Score <- num(1, 7, 0.05, 2)
  d$ASA_CLASS <- codes(1:4, c(0.02, 0.18, 0.6, 0.2))
  d$HTN_YN <- bin(0.88)
  d$HYP_YN <- bin(0.75)
  d$LIPID_YN <- bin(0.9)
  d$PAD_YN <- rep(1, n)
  d$DIAB_YN <- with_na(diabetes, 0.01)
  d$DIAB_TX <- ifelse(diabetes == 1, sample(1:4, n, TRUE), 0)
  d$CAD_YN <- bin(0.4)
  d$ANGINA_YN <- codes(0:2, c(0.85, 0.1, 0.05))
  d$PREV_MI <- bin(0.25)
  d$PREV_PCI <- bin(0.25)
  d$PREV_CABG <- bin(0.15)
  d$COPD_YN <- with_na(copd, 0.01)
  d$TIA_YN <- bin(0.08)
  d$STROKE_YN <- bin(0.1)
  d$CHF_NYHA <- with_na(chf, 0.01)
  d$CKDGRADE <- with_na(ckd, 0.01)
  d$SMOKE_HX <- codes(1:4, c(0.15, 0.45, p$p_current_smoker, 0.04))
  d$SMOKE_PYR <- num(0, 120, 0.27, 0)
  d$SMOK_TX_YN <- codes(0:2)
  prev_infr <- bin(0.4)
  d$PREV_INFR_IND.x <- ifelse(open | both_sides, prev_infr, NA)
  d$PREV_INFR_IND.y <- ifelse(!open | both_sides, prev_infr, NA)
  d$PREV_IND <- bin(0.3)
  d$PREV_INF_IND <- bin(0.5, 0.86)
  d$STATIN <- bin(p$p_statin)
  d$ANTIPLAT <- bin(0.85)
  d$CLOPIDOGREL <- bin(p$p_clopidogrel)
  d$WARFARIN <- bin(0.08)
  d$DOAC <- bin(0.07)
  d$ASPIRIN <- bin(p$p_aspirin)
  d$TICAGRELOR <- bin(0.03)
  d$PRASURGEL <- bin(0.01)
  d$DAPT <- bin(0.2)
  d$CLAUDIC <- bin(0.3, 0.03)
  d$REST_PAIN <- bin(0.5, 0.03)
  d$TISS_LOSS <- bin(0.75, 0.03)
  d$RUTH_4 <- bin(0.25, 0.84)
  d$RUTH_5_6 <- bin(0.75, 0.84)
  d$WAGNER <- codes(0:5, na = 0.03)
  d$DURAT_IND <- codes(1:3, na = 0.03)
  d$CONTRA_LIMB <- bin(0.2)
  d$LESION_LOC <- codes(1:10, na = 0.03)
  d$ALBUMIN.x <- ifelse(open | both_sides, albumin_obs, NA)
  d$ALBUMIN.y <- ifelse(!open | both_sides, albumin_obs, NA)
  d$ALBUMIN_YN <- as.integer(!is.na(albumin_obs))
  d$HEMOGLOBIN <- with_na(ifelse(ifcc, round((hba1c - 2.152) / 0.09148), hba1c), 0.49)
  d$HEMOG_YN <- as.integer(!is.na(d$HEMOGLOBIN))
  d$TOTAL_CHOL <- num(80, 280, 0.67, 0)
  d$LDL <- num(30, 200, 0.68, 0)
  d$LDL[sample(n, 2)] <- 0
  d$HDL <- num(20, 90, 0.67, 0)
  d$TRIGLYCERIDES <- num(50, 400, 0.68, 0)
  d$TRIGLYCERIDES[sample(n, 2)] <- 0
  d$CREAT_YN <- bin(0.98)
  d$CREAT_RES <- with_na(ifelse(creat_umol, round(creat * 88.42), creat), 0.02)
  d$CREAT_UN <- with_na(ifelse(creat_umol, 2, 1), 0.02)
  d$eGFR <- with_na(round(egfr, 1), 0.02)
  d$WIFI.x <- ifelse(open | both_sides, wifi, NA)
  d$WIFI.y <- ifelse(!open | both_sides, wifi, NA)
  d$ISCHEMIA <- with_na(ischemia, 0.06)
  d$INFECTION <- with_na(infection, 0.03)
  d$WIFISTAGE <- codes(1:4, na = 0.04)
  d$VITALS_YN <- bin(0.97)
  d$HEART_R8 <- num(50, 110, 0.03, 0)
  d$SYSTOLIC <- num(90, 190, 0.03, 0)
  d$DIASTOLIC <- num(45, 100, 0.03, 0)
  d$BP_ARM <- codes(1:2)
  d$R_BRACH_PR <- codes(1:2, na = 0.02)
  d$L_BRACH_PR <- codes(1:2, na = 0.02)
  d$SBP_R <- num(90, 190, 0.5, 0)
  d$SBP_L <- num(90, 190, 0.5, 0)
  incompressible <- runif(n) < 0.08
  abi <- ifelse(incompressible, runif(n, 1.41, 2), clip(rnorm(n, p$abi_mean, 0.2), 0, 1.4))
  d$ABIMAX <- with_na(round(abi, 2), 0.32)
  d$TBI <- num(0, 1, 0.59, 2)
  for (v in c("PRESS_DP_IND", "PRESS_PT_IND", "PRESS_DPPT_IND", "PRESS_TOE_IND", "PRESS_DP_NIL", "PRESS_PT_NIL")) {
    d[[v]] <- num(0, 200, 0.4, 0)
  }
  for (v in c("DOPPLER_IND", "DOPPLER_NIL")) d[[v]] <- codes(0:2, na = 0.02)
  for (v in c(outer(c("FEM_P", "POP_P", "POST_P", "DORS_PED"), c("IND", "NIL"), paste, sep = "_"))) d[[v]] <- codes(0:3)
  for (v in c("DP_IND", "PT_IND", "DP_NIL", "PT_NIL", "PVR_IND", "PVR_NIL")) d[[v]] <- codes(1:3, na = 0.02)
  for (v in c("CHEST_IND", "CHEST_NIL", "TMA_IND", "TMA_NIL")) d[[v]] <- num(20, 80, 0.96, 0)
  for (v in c(outer(c("FEM_DOP", "DP_DOP", "PT_DOP"), c("IND", "NIL"), paste, sep = "_"))) {
    d[[v]] <- codes(1:4, na = 0.4)
  }
  for (v in c(outer(c("T1", "T2", "T3", "T4", "T5", "FF", "HF", "AN", "AABK"), c("ULCR", "GANG"), paste, sep = "_"))) {
    d[[v]] <- codes(0:2, c(0.7, 0.2, 0.1), na = 0.25)
  }
  d$TCPO2_YN_IND <- bin(0.1)

  only <- function(x, rows) ifelse(rows, x, NA)
  for (v in c("GREAT_SAPH", "AUTOG_YN", "PROSTH_COND_YN", "NON_CONTIG", "DIST_VEIN", "END_END", "TECH_SUC")) {
    d[[v]] <- only(bin(0.6, 0), open)
  }
  d$DIAM_VEIN <- only(num(2, 6, 0.1), open)
  d$MAX_DIAM <- num(2, 8, 0.14)
  for (v in c("VEIN_HARV", "COND_MANAG")) d[[v]] <- only(codes(1:3, na = 0.25), open)
  d$VEIN_SIDE <- codes(1:2, na = 0.14)
  d$AUTOGEN_YN <- only(bin(0.5, 0.5), open)
  d$AUTOG_VEIN <- only(codes(1:4, na = 0.85), open)
  d$GRAFT_DONOR <- only(codes(1:12, na = 0.02), open)
  d$GRAFT_RECP <- only(codes(1:14), open)
  bypass <- c(
    "BYPASS_CMPSTVN", "AKFMRLPPLTL_BYPS", "BKFMRLPPLTL_BYPS", "FMRLPPLTL_BYPS", "FMRLPDLTBL_BYPS", "PPLTLPDLTBL_BYPS"
  )
  for (v in bypass) {
    d[[v]] <- with_na(ifelse(open, rbinom(n, 1, 0.3), 0), 0.03)
  }
  d$COMP_STUDY <- only(codes(1:4), open)
  # Open-surgery (.x) and endovascular (.y) forms disagree where both were filled in, so both are kept
  d$TRANSFUSION.x <- only(codes(0:4, c(0.7, 0.15, 0.1, 0.03, 0.02)), open | both_sides)
  d$DAYS_ICU.x <- only(num(0, 10, 0, 0), open | both_sides)
  d$BLODD_LOSS <- only(num(50, 1500, 0.02, 0), open)
  d$SEGMENTS <- only(codes(1:6), !open)
  d$SEG_1 <- only(codes(1:18, na = 0.1), !open)
  d$INIT_STENOSIS_1 <- only(num(50, 100, 0.1, 0), !open)
  d$FINAL_STENOSIS_1 <- only(num(0, 50, 0.1, 0), !open)
  d$CTO_1_TREAT <- only(num(0, 300, 0.1, 0), !open)
  d$OSTIAL_1 <- only(bin(0.1, 0.1), !open)
  d$DEVICES_1 <- only(codes(0:6, na = 0.1), !open)
  d$ENDO_TECH_1 <- ifelse(!open, sample(c("PTA", "PTA, stent", "atherectomy", "DCB", ""), n, TRUE), NA)
  d$EMBOLIC_1 <- only(bin(0.05, 0.1), !open)
  for (v in c(
    "SFA_ENDO", "BTK_ENDO", "POPLTL_ENDO", "TBLSPDLS_ENDO", "Artherectomy", "AngioplastyAlone",
    "DrugCtdBlnAngioplasty", "Stents", "DrugElutingstents", "StentGrafts",
    c(outer(c("SFA", "POP", "TBPD"), c("ATH", "PBAA", "PBA", "DCB", "BMS", "DES", "SG"), paste, sep = "_")),
    "ENDO_AORT_THRPY", "IPSIL_CDFA_ENDART", "IPSIL_CDFA_ENDO"
  )) {
    d[[v]] <- with_na(ifelse(open, 0, rbinom(n, 1, 0.3)), 0.03)
  }
  d$SUCCESS_1 <- only(bin(0.85, 0), !open)
  d$TRANSFUSION.y <- only(codes(0:4, c(0.9, 0.05, 0.03, 0.01, 0.01)), !open | both_sides)
  d$DAYS_ICU.y <- only(num(0, 5, 0, 0), !open | both_sides)
  d$CONTRAST <- only(num(20, 300, 0.05, 0), !open)
  d$CO2_CONTRAST <- only(bin(0.1, 0), !open)
  d$SHEATH <- only(codes(4:8), !open)
  d$REVASC_PROC <- only(bin(0.15, 0), !open)
  d$PROC_FRMRND <- codes(1:3)
  d$CONCUR_AMP <- bin(0.05)
  d$REINT_ADJ_YN <- bin(0.2)
  d$WHO_DECIDE <- codes(1:3, na = 0.95)
  d$UNUSED_FIELD <- rep(NA, n)
  d$EVNT_AAAIOS <- amp_event
  d$TMTO_AAAIOS <- amp_time
  d$EVNT_DTH <- death_event
  d$TMTO_DTHYR <- death_time
  d$EVNT_REINTIOS <- reint_event
  d$TMTO_REINTIOS <- reint_time
  d$EVNT_MALEDTH <- male_event
  d$TMTO_MALEDTH <- male_time

  # Near-empty forms: four or more comorbidity boxes blank
  empty <- sample(n, 3)
  d[empty, c("HTN_YN", "DIAB_YN", "CAD_YN", "COPD_YN", "PREV_MI")] <- NA

  randomization <- data.frame(
    Random_Subnum = d$Random_Subnum, VISNAME = "Screening", PAGENAME = "Randomization",
    VISITID = 20, VISITSEQ = 0, PAGESEQ = 0,
    FACTOR_A = sample(1:2, n, TRUE), FACTOR_B = 1 + infrapopliteal, INDEX_LEG = sample(1:2, n, TRUE)
  )
  # As in the trial export, two patients treated in both arms have a second procedure record
  twice <- head(setdiff(seq_len(n), empty), 2)
  d <- rbind(d, d[twice, ])
  list(export = d, randomization = randomization)
}
