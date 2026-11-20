#' Pages of the BEST-CLI SAS datasets that make up the export, in merge order.
#' `baseline_visit` keeps the Baseline row only; `keep` and `drop` select columns.
#' Not used: med_hist (its processed form is in baseline), reint, vitalsda, proc_test.
export_pages <- function() {
  visit_cols <- c("VISNAME", "PAGENAME", "VISITID", "VISITSEQ", "PAGESEQ")
  outcomes <- c(
    "TMTO_AAAIOS", "EVNT_AAAIOS", "TMTO_REINTIOS", "EVNT_REINTIOS",
    "TMTO_MALEDTH", "EVNT_MALEDTH", "EVNT_DTH", "TMTO_DTHYR"
  )
  # Visit bookkeeping stays only on the two procedure pages, where it tells their records apart
  list(
    baseline = list(),
    timetoevent_outcomes = list(keep = outcomes),
    labs = list(baseline_visit = TRUE, drop = visit_cols),
    vasc_exam = list(baseline_visit = TRUE, drop = c(visit_cols, "ASSESS_FRMRND")),
    vitals = list(baseline_visit = TRUE, drop = visit_cols),
    hem_assess = list(baseline_visit = TRUE, drop = c(visit_cols, "ASSESS_FRMRND")),
    open_surg = list(),
    endo = list(),
    vein_map = list(drop = c(visit_cols, "EXAM_FRMRND")),
    vasc_hx = list(drop = c(visit_cols, "ASSESS_FRMRND")),
    index_sum = list(drop = "PROC_FRMRND"),
    proc_indexvst = list()
  )
}

#' Read one SAS page and apply its row and column selection.
read_export_page <- function(sas_dir, name, spec) {
  path <- file.path(sas_dir, paste0(name, ".sas7bdat"))
  if (!file.exists(path)) stop("Missing ", path)
  df <- as.data.frame(haven::zap_label(haven::zap_labels(haven::zap_formats(haven::read_sas(path)))))
  if (isTRUE(spec$baseline_visit)) {
    df <- df[df$VISNAME == "Baseline", , drop = FALSE]
    if (anyDuplicated(df$Random_Subnum)) stop(name, ": more than one Baseline row for some patients")
  }
  if (!is.null(spec$keep)) {
    miss <- setdiff(spec$keep, names(df))
    if (length(miss)) stop(name, ": missing columns ", paste(miss, collapse = ", "))
    df <- df[, c("Random_Subnum", spec$keep), drop = FALSE]
  }
  if (!is.null(spec$drop)) df <- df[, setdiff(names(df), spec$drop), drop = FALSE]
  # merge() re-sorts rows, so a patient's second record on a page keeps its SAS order through this column
  if (anyDuplicated(df$Random_Subnum)) {
    df[[paste0(".seq_", name)]] <- ave(seq_len(nrow(df)), df$Random_Subnum, FUN = seq_along)
  }
  df
}

#' Assemble FullAugmentedData and Randomization from the SAS datasets.
#' Every page is left-joined onto the baseline summary; column order and the .x (open) / .y (endo)
#' suffixes follow from the merge order in `export_pages()`.
#' @return List: `data` (one row per procedure record) and `randomization` (one row per patient).
build_trial_export <- function(sas_dir) {
  specs <- export_pages()
  pages <- Map(function(name, spec) read_export_page(sas_dir, name, spec), names(specs), specs)
  df <- Reduce(function(x, y) merge(x, y, by = "Random_Subnum", all.x = TRUE), pages)
  seq_cols <- grep("^\\.seq_", names(df), value = TRUE)
  df <- df[do.call(order, unname(as.list(df[, c("Random_Subnum", seq_cols), drop = FALSE]))), , drop = FALSE]
  df <- df[, setdiff(names(df), seq_cols), drop = FALSE]
  rownames(df) <- NULL
  rand_path <- file.path(sas_dir, "rand.sas7bdat")
  if (!file.exists(rand_path)) stop("Missing ", rand_path)
  list(data = df, randomization = haven::read_sas(rand_path))
}
