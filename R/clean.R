#' Round to the 10 significant digits that Excel's General format keeps, halves away from zero.
#' The original CLI_clean.csv passed through Excel; matching its rounding reproduces the analysed data.
round_like_excel <- function(x) {
  k <- !is.na(x) & x != 0
  d <- 10 - pmax(1, floor(log10(abs(x[k]))) + 1)
  # signif(, 15) settles binary near-ties as Excel's 15-digit decimal does
  z <- floor(signif(abs(x[k]) * 10^d, 15) + 0.5)
  x[k] <- sign(x[k]) * as.numeric(sprintf("%.*f", d, z / 10^d))
  x
}

#' For a .x/.y column pair: rows where both are observed, and how many of those disagree.
pair_concordance <- function(a, b, tol = 1e-8) {
  both <- !is.na(a) & !is.na(b)
  disagree <- if (is.numeric(a) && is.numeric(b)) {
    sum(abs(a[both] - b[both]) > tol)
  } else {
    sum(as.character(a[both]) != as.character(b[both]))
  }
  list(n_overlap = sum(both), n_disagree = disagree)
}

#' Columns kept in the cleaned data, by role.
clean_columns <- function() {
  list(
    baseline = c(
      "AGE", "GENDER", "HISPANIC", "BMI", "HEIGHT", "WEIGHT",
      "AMB_STAT3CAT", "LIVE_STAT", "vascu_qol_Score", "ASA_CLASS",
      "HTN_YN", "HYP_YN", "DIAB_YN", "DIAB_TX", "CAD_YN", "ANGINA_YN", "PREV_MI", "PREV_PCI", "PREV_CABG",
      "COPD_YN", "TIA_YN", "STROKE_YN", "CHF_NYHA", "ESKD", "SMOKE_HX", "SMOKE_PYR", "SMOK_TX_YN",
      "PREV_INFR_IND", "PREV_IND", "PREV_INF_IND",
      "STATIN", "ANTIPLAT", "CLOPIDOGREL", "WARFARIN", "DOAC", "ASPIRIN", "TICAGRELOR", "PRASURGEL", "DAPT",
      "CLAUDIC", "REST_PAIN", "TISS_LOSS", "RUTH_4", "RUTH_5_6", "WAGNER", "DURAT_IND", "CONTRA_LIMB", "LESION_LOC",
      "ALBUMIN", "HBA1C", "TOTAL_CHOL", "LDL", "HDL", "TRIGLYCERIDES", "CREAT_RES", "eGFR",
      "Wound", "ISCHEMIA", "INFECTION",
      "VITALS_YN", "HEART_R8", "SYSTOLIC", "DIASTOLIC", "R_BRACH_PR", "L_BRACH_PR", "SBP_R", "SBP_L",
      "ABIMAX", "ABI_INCOMPRESSIBLE", "TBI", "PRESS_DP_IND", "PRESS_PT_IND", "PRESS_DPPT_IND", "PRESS_TOE_IND",
      "DOPPLER_IND", "PRESS_DP_NIL", "DOPPLER_NIL",
      "FEM_P_IND", "FEM_P_NIL", "POP_P_IND", "POP_P_NIL", "POST_P_IND", "POST_P_NIL", "DORS_PED_IND", "DORS_PED_NIL",
      "DP_IND", "PT_IND", "DP_NIL", "PT_NIL", "CHEST_IND", "CHEST_NIL", "TMA_IND", "TMA_NIL", "PVR_IND", "PVR_NIL",
      "FEM_DOP_IND", "FEM_DOP_NIL", "DP_DOP_IND", "DP_DOP_NIL", "PT_DOP_IND", "PT_DOP_NIL",
      paste0(rep(c("T1", "T2", "T3", "T4", "T5", "FF", "HF", "AN", "AABK"), each = 2), "_", c("ULCR", "GANG"))
    ),
    procedural = c(
      "GREAT_SAPH", "AUTOG_YN", "PROSTH_COND_YN", "NON_CONTIG", "DIST_VEIN", "END_END", "DIAM_VEIN", "MAX_DIAM",
      "VEIN_HARV", "COND_MANAG", "VEIN_SIDE", "AUTOGEN_YN", "AUTOG_VEIN", "GRAFT_DONOR", "GRAFT_RECP",
      "BYPASS_CMPSTVN", "AKFMRLPPLTL_BYPS", "BKFMRLPPLTL_BYPS", "FMRLPPLTL_BYPS", "FMRLPDLTBL_BYPS",
      "PPLTLPDLTBL_BYPS", "TECH_SUC", "COMP_STUDY", "TRANSFUSION.x", "DAYS_ICU.x", "BLODD_LOSS",
      "SEGMENTS", "SEG_1", "INIT_STENOSIS_1", "FINAL_STENOSIS_1", "CTO_1_TREAT", "OSTIAL_1", "DEVICES_1",
      "ENDO_TECH_1", "EMBOLIC_1", "SFA_ENDO", "BTK_ENDO", "POPLTL_ENDO", "TBLSPDLS_ENDO",
      "Artherectomy", "AngioplastyAlone", "DrugCtdBlnAngioplasty", "Stents", "DrugElutingstents", "StentGrafts",
      "SUCCESS_1", "TRANSFUSION.y", "DAYS_ICU.y", "CONTRAST", "CO2_CONTRAST", "SHEATH",
      paste0(rep(c("SFA", "POP", "TBPD"), each = 7), "_", c("ATH", "PBAA", "PBA", "DCB", "BMS", "DES", "SG")),
      "REVASC_PROC", "ENDO_AORT_THRPY", "IPSIL_CDFA_ENDART", "IPSIL_CDFA_ENDO"
    ),
    outcome = c(
      "EVNT_AAAIOS", "TMTO_AAAIOS", "EVNT_DTH", "TMTO_DTHYR",
      "EVNT_REINTIOS", "TMTO_REINTIOS", "EVNT_MALEDTH", "TMTO_MALEDTH",
      "cr_status", "cr_time", "amp_event", "amp_time", "death_event", "death_time", "afs_event", "afs_time"
    )
  )
}

#' Clean the raw BEST-CLI export: harmonize units, apply the record exclusions, derive the outcomes.
#' @param df `FullAugmentedData.csv`; `randomization`: `Randomization.csv`.
#' @return List: `data` (one row per patient) and `manifest` (every column and row decision).
clean_trial_data <- function(df, randomization) {
  umol_per_mgdl <- 88.42
  cm_per_inch <- 2.54
  kg_per_lb <- 0.45359237

  log <- new.env()
  log$steps <- list()
  log_step <- function(stage, column, action, detail = "") {
    log$steps[[length(log$steps) + 1]] <- data.frame(stage = stage, column = column, action = action, detail = detail)
  }
  drop_cols <- function(df, cols, stage, detail) {
    cols <- intersect(cols, names(df))
    for (c in cols) log_step(stage, c, "drop", detail)
    df[, !(names(df) %in% cols), drop = FALSE]
  }

  all_na <- names(df)[vapply(df, function(x) all(is.na(x)), logical(1))]
  df <- drop_cols(df, all_na, "1_drop_all_na", "entirely missing")

  # The export merged the open-surgery and endovascular tables, leaving .x/.y pairs.
  # Concordant pairs are combined; discordant ones are open (.x) vs endo (.y) fields and both are kept.
  suffixed <- grep("\\.(x|y)$", names(df), value = TRUE)
  for (b in unique(sub("\\.(x|y)$", "", suffixed))) {
    xc <- paste0(b, ".x")
    yc <- paste0(b, ".y")
    have <- intersect(c(xc, yc), names(df))
    if (length(have) == 1) {
      if (!(b %in% names(df))) {
        names(df)[names(df) == have] <- b
        log_step("1_xy", have, "rename_to_base", paste0("only one side present -> ", b))
      }
    } else {
      cc <- pair_concordance(df[[xc]], df[[yc]])
      if (cc$n_disagree == 0 && !(b %in% names(df))) {
        df[[b]] <- df[[xc]]
        df[[b]][is.na(df[[b]])] <- df[[yc]][is.na(df[[b]])]
        df[[xc]] <- NULL
        df[[yc]] <- NULL
        log_step("1_xy", paste(xc, yc, sep = " + "), "coalesce", sprintf("concordant; overlap=%d", cc$n_overlap))
      } else {
        log_step(
          "1_xy", paste(xc, yc, sep = " + "), "FLAG_discordant",
          sprintf("overlap=%d disagree=%d -> kept both (open .x / endo .y)", cc$n_overlap, cc$n_disagree)
        )
      }
    }
  }

  as_text <- function(x) paste(ifelse(is.na(x), "\001NA\001", as.character(x)), collapse = "\002")
  signature <- vapply(df, as_text, character(1))
  first_match <- match(signature, signature)
  for (i in which(seq_along(signature) != first_match)) {
    log_step("1_exact_dup", names(df)[i], "drop", paste0("byte-identical to ", names(df)[first_match[i]]))
  }
  df <- df[, seq_along(signature) == first_match, drop = FALSE]

  # ISCHEMIA was exported 1-4; the SVS WIfI grade is 0-3
  if (all(df$ISCHEMIA[!is.na(df$ISCHEMIA)] %in% 1:4)) {
    df$ISCHEMIA <- df$ISCHEMIA - 1
    log_step("2_ischemia", "ISCHEMIA", "reindex", "1-4 -> 0-3 (SVS norm)")
  }
  # WIFI holds the wound component only, not the composite stage
  names(df)[names(df) == "WIFI"] <- "Wound"
  log_step("2_rename", "WIFI -> Wound", "rename", "WIfI wound grade 0-3")
  df <- drop_cols(df, "WIFISTAGE", "2_drop", "composite stage; components kept")
  umol <- !is.na(df$CREAT_UN) & df$CREAT_UN == 2
  df$CREAT_RES[umol] <- df$CREAT_RES[umol] / umol_per_mgdl
  df$CREAT_UN[!is.na(df$CREAT_UN)] <- 1
  log_step("2_creat", "CREAT_RES", "unit_convert", sprintf("%d umol/L -> mg/dL", sum(umol)))
  df <- drop_cols(df, grep("_FRMRND", names(df), value = TRUE), "2_frmrnd", "measurement timing (form round)")

  df$eGFR <- pmin(pmax(df$eGFR, 5), 120)
  log_step("3_egfr", "eGFR", "winsorize", "clipped to [5, 120]")
  df$ESKD <- ifelse(is.na(df$CKDGRADE), NA_integer_, as.integer(df$CKDGRADE == 5))
  log_step("3_renal", "ESKD", "derive", "= CKDGRADE == 5 (kidney transplant; see docs/decisions.md)")
  df <- drop_cols(df, "CKDGRADE", "3_renal", "re-binning of eGFR")

  df <- drop_cols(df, c("WHO_DECIDE", "ASSIGN_REAS"), "4_drop", "sparse crossover-decision metadata")
  df <- drop_cols(df, "RAND_ASSIGN", "5_trt", "duplicate of ARM")
  df <- drop_cols(df, "ACTTRT", "5_trt", "superseded by ASSIGNED")
  df <- drop_cols(df, "BP_ARM", "5_trt", "side of blood-pressure measurement")

  comorbidities <- intersect(c(
    "HTN_YN", "HYP_YN", "DIAB_YN", "CAD_YN", "COPD_YN", "TIA_YN", "STROKE_YN", "LIPID_YN",
    "PREV_MI", "PREV_PCI", "PREV_CABG", "STATIN", "ASPIRIN", "CLOPIDOGREL", "WARFARIN", "DOAC",
    "TICAGRELOR", "PRASURGEL", "ESKD", "CHF_NYHA", "PREV_INFR_IND", "PREV_IND", "PAD_YN"
  ), names(df))
  # Records missing 4 or more comorbidity flags are near-empty forms
  near_empty <- which(rowSums(is.na(df[, comorbidities, drop = FALSE])) >= 4)
  for (id in df$Random_Subnum[near_empty]) log_step("6_exclude_rows", as.character(id), "drop_row", "near-empty record")
  if (length(near_empty)) df <- df[-near_empty, , drop = FALSE]
  # An unchecked comorbidity box is recorded as missing; treat it as absent
  for (c in setdiff(comorbidities, "PAD_YN")) {
    k <- sum(is.na(df[[c]]))
    if (k > 0) {
      df[[c]][is.na(df[[c]])] <- 0
      log_step("6_na_to_zero", c, "recode_NA_0", sprintf("%d NA -> 0", k))
    }
  }
  df <- drop_cols(df, c("REINT_ADJ_YN", "CREAT_UN"), "6_drop", "not a comorbidity flag")
  df$HEIGHT[!is.na(df$HEIGHT) & df$HEIGHT == 0] <- NA
  inch <- !is.na(df$HEIGHT_UN) & df$HEIGHT_UN == 2
  df$HEIGHT[inch] <- df$HEIGHT[inch] * cm_per_inch
  df$HEIGHT_UN <- NULL
  log_step("6_height", "HEIGHT", "unit_convert", "inches -> cm; 0 -> NA")
  lb <- !is.na(df$WEIGHT_UN) & df$WEIGHT_UN == 2
  df$WEIGHT[lb] <- df$WEIGHT[lb] * kg_per_lb
  df$WEIGHT_UN <- NULL
  log_step("6_weight", "WEIGHT", "unit_convert", "lb -> kg")
  message(sprintf("Excluded %d near-empty records", length(near_empty)))

  df <- drop_cols(df, "PAD_YN", "7_drop", "uninformative in an all-PAD cohort")
  df <- drop_cols(
    df, c("CREAT_YN", "ALBUMIN_YN", "HEMOG_YN", "LIPID_YN"), "8_drop_lab_YN", "test-performed flag; value kept"
  )
  df <- drop_cols(df, "AMB_STAT", "8_drop", "4-level version; AMB_STAT3CAT kept")

  # HEMOGLOBIN holds HbA1c, partly in IFCC mmol/mol; values above 20 are converted to NGSP %
  names(df)[names(df) == "HEMOGLOBIN"] <- "HBA1C"
  ifcc <- !is.na(df$HBA1C) & df$HBA1C > 20
  df$HBA1C[ifcc] <- 0.09148 * df$HBA1C[ifcc] + 2.152
  log_step("9_hba1c", "HEMOGLOBIN -> HBA1C", "rename+harmonize", sprintf("IFCC -> NGSP for %d", sum(ifcc)))
  for (c in c("LDL", "TRIGLYCERIDES")) {
    zero <- !is.na(df[[c]]) & df[[c]] == 0
    if (any(zero)) {
      df[[c]][zero] <- NA
      log_step("9_lipids", c, "zero_to_NA", sprintf("%d impossible zeros -> NA", sum(zero)))
    }
  }

  df$ABI_INCOMPRESSIBLE <- ifelse(is.na(df$ABIMAX), NA_integer_, as.integer(df$ABIMAX > 1.4))
  log_step("10_abi", "ABI_INCOMPRESSIBLE", "derive", "= ABIMAX > 1.4")
  df <- drop_cols(
    df, c("TCPO2_YN_IND", "TCPO2_YN_NIL", "PRESS_PT_NIL", "PRESS_TOE_NIL", "TOE_NIL", "TOE_IND"),
    "10_drop_hemo", "no-value indicator or redundant contralateral measure"
  )
  df <- drop_cols(df, "CONCUR_AMP", "11_drop", "recorded after baseline; leaks the amputation outcome")

  # The export holds one record per index procedure, so a patient treated in both arms appears twice
  two_treatments <- sort(unique(df$Random_Subnum[duplicated(df$Random_Subnum)]))
  df <- df[!(df$Random_Subnum %in% two_treatments), ]
  for (id in two_treatments) log_step("12_drop_row", as.character(id), "drop_row", "received two treatments")

  amp <- df$EVNT_AAAIOS == 1
  death <- df$EVNT_DTH == 1
  df$cr_status <- ifelse(amp, 1L, ifelse(death, 2L, 0L))
  df$cr_time <- ifelse(amp, df$TMTO_AAAIOS, ifelse(death, df$TMTO_DTHYR, pmin(df$TMTO_AAAIOS, df$TMTO_DTHYR)))
  df$amp_event <- as.integer(df$EVNT_AAAIOS)
  df$amp_time <- df$TMTO_AAAIOS
  df$death_event <- as.integer(df$EVNT_DTH)
  df$death_time <- df$TMTO_DTHYR
  df$afs_event <- as.integer(amp | death)
  df$afs_time <- df$cr_time
  log_step("12_targets", "cr_status/cr_time/amp/death/afs", "derive", "competing risks and amputation-free survival")

  # One patient's amputation postdates their vital-status censoring; keep the
  # amputation and extend death censoring to it
  late <- which(df$amp_event == 1 & df$death_event == 0 & df$death_time < df$amp_time)
  df$TMTO_DTHYR[late] <- df$amp_time[late]
  df$death_time[late] <- df$amp_time[late]
  for (id in df$Random_Subnum[late]) {
    log_step("14_anomaly", as.character(id), "fix_censor", "death censoring extended to amputation time")
  }

  from_sheet <- c("INDEX_LEG", "FACTOR_A", "FACTOR_B")
  df <- merge(df, randomization[, c("Random_Subnum", from_sheet)], by = "Random_Subnum", all.x = TRUE, sort = FALSE)
  for (c in from_sheet) log_step("15_randomization", c, "merge", "joined from Randomization.csv")

  cols <- clean_columns()
  randomization_cols <- c(intersect(c("ARM", "ASSIGNED", "STRATUM", "COHORT", "SSGSV_YN"), names(df)), from_sheet)
  missing_cols <- setdiff(c(cols$baseline, cols$procedural), names(df))
  if (length(missing_cols)) stop("Expected columns absent: ", paste(missing_cols, collapse = ", "))
  keep <- c("Random_Subnum", randomization_cols, cols$baseline, cols$procedural, cols$outcome)
  df <- df[, keep]
  for (v in names(df)[vapply(df, is.double, logical(1))]) df[[v]] <- round_like_excel(df[[v]])
  log_step("16_round", "all numeric", "round", "10 significant digits, as Excel stored the original CLI_clean.csv")

  list(data = df, manifest = do.call(rbind, log$steps))
}
