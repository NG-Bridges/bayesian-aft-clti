#' Look up one pooled coefficient row; anything other than exactly one row is an error.
coef_row <- function(co, model, cause, variable) {
  z <- co[co$model == model & co$cause == cause & co$variable == variable, ]
  if (nrow(z) != 1) stop(sprintf("expected one row for %s / %s / %s, found %d", model, cause, variable, nrow(z)))
  z
}

coef_cell <- function(co, model, cause, variable, d = 2, bold_ok = TRUE) {
  z <- coef_row(co, model, cause, variable)
  tex_cri(z$post_median, z$lo95, z$hi95, d, bold = bold_ok && excludes_zero(z$lo95, z$hi95))
}

pooling_note <- paste(
  "Posterior median (95\\% equal-tailed credible interval) of the time-ratio coefficient on the log-time scale;",
  "positive values lengthen time to the event. Draws were pooled by stacking 1,000 posterior draws from each of",
  "the 50 imputed-dataset fits."
)
design_note <- "Design indicators: Cohort 2 = no adequate single-segment great saphenous vein (GSV) conduit; endovascular arm vs.\\ surgical bypass (randomized assignment)."

#' WAIC/DIC3 table values, with a combined SE for each difference from the single-component model.
#' @return Data frame: one row per cause and model.
ic_values <- function(ics, pdiff, m) {
  rows <- list()
  for (cs in names(cause_labels())) {
    for (mdl in names(model_labels())) {
      s <- ics[ics$model == mdl & ics$cause == cs, ]
      if (mdl == "AFT_single") {
        d <- data.frame(dW = NA, se_w = NA, sd_a = NA, dD = NA)
      } else {
        # Paired differences are (first - second); re-signed so negative favours the model
        p <- pdiff[pdiff$cause == cs & pdiff$comparison == paste("AFT_single -", mdl), ]
        if (nrow(p) != 1) stop("paired row missing for ", mdl, " / ", cs)
        d <- data.frame(dW = -p$d_waic, se_w = p$d_waic_se_within, sd_a = p$d_waic_sd_across, dD = -p$d_dic3)
      }
      se_c <- sqrt(d$se_w^2 + (1 + 1 / m) * d$sd_a^2)
      rows[[length(rows) + 1]] <- data.frame(
        cause = cs, model = mdl, waic = s$waic, waic_sd = s$waic_sd,
        p_waic = s$p_waic, dWAIC = d$dW, se_within = d$se_w, sd_across = d$sd_a, se_combined = se_c,
        ratio = d$dW / se_c, dDIC3 = d$dD, n_flagged = s$n_flagged_obs, max_p_waic_i = s$max_p_waic_i
      )
    }
  }
  do.call(rbind, rows)
}

table_ic_tex <- function(ic, m) {
  body <- character(0)
  for (cs in names(cause_labels())) {
    body <- c(body, if (length(body)) "\\addlinespace", tex_panel_head(cause_labels()[cs], 7))
    for (r in which(ic$cause == cs)) {
      z <- ic[r, ]
      ref <- is.na(z$dWAIC)
      body <- c(body, sprintf(
        "\\quad %s & %s (%s) & %s & %s & %s & %s & %s \\\\",
        model_labels()[z$model], tex_num(z$waic, 1), tex_num(z$waic_sd, 1), tex_num(z$p_waic, 1),
        if (ref) "ref." else tex_num(z$dWAIC, 1), if (ref) "---" else tex_num(z$se_combined, 2),
        if (ref) "---" else tex_num(z$ratio, 1), if (ref) "ref." else tex_num(z$dDIC3, 1)
      ))
    }
  }
  c(
    tex_header(),
    "\\begin{table}[htbp]", "\\centering", "\\footnotesize",
    "\\begin{threeparttable}",
    "\\caption{Information criteria across the 50 imputed datasets.}\\label{tab:ic}",
    "\\begin{tabular}{@{}lcccccc@{}}", "\\toprule",
    " & WAIC, mean (SD)\\tnote{a} & $p_{\\mathrm{WAIC}}$ & $\\Delta$WAIC\\tnote{b} & SE$_{\\Delta}$\\tnote{c} & $\\Delta$WAIC/SE$_{\\Delta}$ & $\\Delta$DIC$_3$\\tnote{b} \\\\",
    "\\midrule", body, "\\bottomrule", "\\end{tabular}",
    "\\begin{tablenotes}[flushleft]\\scriptsize",
    "\\item[a] Mean and standard deviation across the 50 imputed datasets; lower is better.",
    "\\item[b] Model minus single-component, computed within each imputed dataset and averaged; negative values favor the model.",
    sprintf("\\item[c] Combined standard error, $\\sqrt{\\mathrm{SE}_{\\mathrm{within}}^2 + (1 + 1/m)\\,\\mathrm{SD}_{\\mathrm{across}}^2}$ with $m = %d$: SE$_{\\mathrm{within}}$ is the patient-level standard error of the paired difference, averaged over datasets, and SD$_{\\mathrm{across}}$ is the standard deviation of the difference across datasets. $|\\Delta\\mathrm{WAIC}/\\mathrm{SE}_{\\Delta}| < 2$ indicates models that are not distinguishable.", m),
    "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
  )
}

#' Mixture coefficients by class (relabelled), with the paired Class 2 - Class 1 difference.
table_coef_mixture_tex <- function(co, mxd) {
  design <- design_labels()
  clin <- clinical_labels()
  diff_cell <- function(cs, v, d = 2) {
    z <- mxd[mxd$cause == cs & mxd$variable == v, ]
    if (nrow(z) != 1) {
      return("---")
    }
    tex_cri(z$median, z$lo95, z$hi95, d, bold = excludes_zero(z$lo95, z$hi95))
  }
  cells <- function(cs, v, d = 2, bold_ok = TRUE) {
    c(
      coef_cell(co, "AFT_mixture", cs, paste0(v, "_c1"), d, bold_ok),
      coef_cell(co, "AFT_mixture", cs, paste0(v, "_c2"), d, bold_ok),
      diff_cell(cs, v, d)
    )
  }
  body <- character(0)
  for (cs in names(cause_labels())) {
    body <- c(
      body, if (length(body)) "\\addlinespace",
      tex_panel_head(sprintf("%s. %s", LETTERS[match(cs, names(cause_labels()))], cause_labels()[cs]), 4),
      tex_group_head("Design variables (0/1; weakly informative prior)", 4)
    )
    for (v in names(design)) body <- c(body, tex_row(design[v], cells(cs, v)))
    body <- c(body, tex_group_head("Clinical features (per SD; regularized horseshoe prior)", 4))
    for (v in names(clin)) body <- c(body, tex_row(clin[v], cells(cs, v)))
    body <- c(body, tex_group_head("Hazard", 4), tex_row("Weibull shape, $k$", cells(cs, "k_shape", 3, FALSE)))
  }
  notes <- paste(
    pooling_note,
    "Component labels were resolved within each fit by Stephens' algorithm and aligned",
    "across fits by ordering components on their allocation-weighted mean eGFR, so Class 1 is the lower-eGFR component",
    "in every dataset. Clinical features were standardized; binary clinical features are therefore per SD of the",
    "indicator. Bold: interval excludes zero.",
    "The difference column is computed draw by draw within each fit and then stacked.",
    design_note
  )
  tex_longtable(
    3, " & Class 1 & Class 2 & Class 2 \\textminus{} Class 1 \\\\",
    "Gated mixture of Weibull AFT models: component-specific coefficients and hazard shapes.",
    "tab:coef-mixture", body, notes
  )
}

table_coef_single_frailty_tex <- function(co) {
  design <- design_labels()
  clin <- clinical_labels()
  cells <- function(cs, v) c(coef_cell(co, "AFT_single", cs, v), coef_cell(co, "AFT_frailty", cs, v))
  body <- character(0)
  for (cs in names(cause_labels())) {
    body <- c(
      body, if (length(body)) "\\addlinespace",
      tex_panel_head(sprintf("%s. %s", LETTERS[match(cs, names(cause_labels()))], cause_labels()[cs]), 3),
      tex_group_head("Design variables (0/1; weakly informative prior)", 3)
    )
    for (v in names(design)) body <- c(body, tex_row(design[v], cells(cs, v)))
    body <- c(body, tex_group_head("Clinical features (per SD; regularized horseshoe prior)", 3))
    for (v in names(clin)) body <- c(body, tex_row(clin[v], cells(cs, v)))
    body <- c(
      body, tex_group_head("Hazard", 3),
      tex_row("Weibull shape, $k$", c(
        coef_cell(co, "AFT_single", cs, "k_shape", 3, FALSE),
        coef_cell(co, "AFT_frailty", cs, "k_shape", 3, FALSE)
      )),
      tex_row("Frailty SD, $\\sigma_u$", c("---", coef_cell(co, "AFT_frailty", cs, "sigma_u", 3, FALSE)))
    )
  }
  notes <- paste(
    pooling_note,
    "The frailty model integrates a log-normal subject frailty out of the likelihood",
    "(25-node Gauss--Hermite quadrature), so its coefficients are population-averaged. With one event time per patient,",
    "$\\sigma_u$ and $k$ are only weakly separately identified. Clinical features were standardized; binary clinical",
    "features are therefore per SD of the indicator. Bold: interval excludes zero.",
    design_note
  )
  tex_longtable(
    2, " & Single-component & Frailty \\\\",
    "Single-component and frailty Weibull AFT models: coefficients and hazard parameters.",
    "tab:coef-single-frailty", body, notes
  )
}

#' Convergence summary per model and cause; the mixture gets label-invariant and relabelled columns.
convergence_values <- function(dg, cfg) {
  dg$converged <- as.logical(dg$converged)
  conv_col <- function(z, which) {
    if (which == "lp") {
      rh <- z$rhat_lp
      es <- z$ess_lp
    } else {
      rh <- z$max_rhat
      es <- z$min_ess
    }
    pass <- rh < cfg$rhat_max & es > cfg$ess_min & z$div_frac <= cfg$div_frac_max & z$treedepth_hits == 0
    c(
      n = nrow(z), pass = sum(pass), rhat_max = max(rh), rhat_med = median(rh),
      n101 = sum(rh >= 1.01), n102 = sum(rh >= 1.02), ess_min = min(es), ess_med = median(es),
      n_div = sum(z$n_divergent > 0), div_max = 100 * max(z$div_frac), td = sum(z$treedepth_hits)
    )
  }
  cv <- list()
  for (cs in names(cause_labels())) {
    s <- dg[dg$model == "AFT_single" & dg$cause == cs, ]
    f <- dg[dg$model == "AFT_frailty" & dg$cause == cs, ]
    x <- dg[dg$model == "AFT_mixture" & dg$cause == cs, ]
    cv[[cs]] <- cbind(SC = conv_col(s, "par"), FM = conv_col(f, "par"), MM_lp = conv_col(x, "lp"), MM_rel = conv_col(x, "par"))
    # The recomputed pass counts must equal the converged flags written by 02_fit_models.R
    stopifnot(
      cv[[cs]]["pass", "SC"] == sum(s$converged), cv[[cs]]["pass", "FM"] == sum(f$converged),
      cv[[cs]]["pass", "MM_rel"] == sum(x$converged)
    )
  }
  cv
}

table_convergence_tex <- function(cv, cfg) {
  fmt_row <- function(lab, v, f, span_sampler = FALSE) {
    cells <- vapply(v, f, character(1))
    if (span_sampler) cells <- c(cells[1:2], sprintf("\\multicolumn{2}{c}{%s}", cells[3]))
    sprintf("\\quad %s & %s \\\\", lab, paste(cells, collapse = " & "))
  }
  int <- function(x) as.character(as.integer(x))
  f3 <- function(x) formatC(x, format = "f", digits = 3)
  big <- function(x) formatC(round(x), format = "d", big.mark = ",")
  body <- character(0)
  for (cs in names(cause_labels())) {
    m <- cv[[cs]]
    n <- m["n", 1]
    body <- c(
      body, if (length(body)) "\\addlinespace", tex_panel_head(cause_labels()[cs], 5),
      fmt_row("Fits meeting all criteria\\tnote{a}", m["pass", ], function(x) sprintf("%d/%d", as.integer(x), as.integer(n))),
      fmt_row("Largest $\\hat{R}$ (worst fit)", m["rhat_max", ], f3),
      fmt_row("Median of per-fit largest $\\hat{R}$", m["rhat_med", ], f3),
      fmt_row("Fits with $\\hat{R} \\geq 1.01$", m["n101", ], int),
      fmt_row("Fits with $\\hat{R} \\geq 1.02$", m["n102", ], int),
      fmt_row("Smallest ESS (worst fit)", m["ess_min", ], big),
      fmt_row("Median of per-fit smallest ESS", m["ess_med", ], big),
      fmt_row("Fits with $\\geq$1 divergent transition\\tnote{b}", m["n_div", 1:3], int, TRUE),
      fmt_row("Largest share of divergent transitions (\\%)\\tnote{b}", m["div_max", 1:3], function(x) formatC(x, format = "f", digits = 2), TRUE),
      fmt_row("Fits reaching maximum tree depth\\tnote{b}", m["td", 1:3], int, TRUE)
    )
  }
  c(
    tex_header(),
    "\\begin{table}[htbp]", "\\centering", "\\footnotesize",
    "\\begin{threeparttable}",
    "\\caption{MCMC convergence across the 50 imputed-dataset fits per model and cause.}\\label{tab:convergence}",
    "\\begin{tabular}{@{}lcccc@{}}", "\\toprule",
    " & & & \\multicolumn{2}{c}{Gated mixture} \\\\", "\\cmidrule(l){4-5}",
    " & Single-component & Frailty & Label-invariant\\tnote{c} & Relabeled\\tnote{d} \\\\",
    "\\midrule", body, "\\bottomrule", "\\end{tabular}",
    "\\begin{tablenotes}[flushleft]\\scriptsize",
    sprintf("\\item[a] $\\hat{R} < %.2f$, ESS $> %d$, divergent transitions $\\leq %.1f\\%%$ of draws, and no draws at the maximum tree depth. Each fit used 4 chains of 1,500 warmup and 1,500 retained iterations.", cfg$rhat_max, cfg$ess_min, 100 * cfg$div_frac_max),
    "\\item[b] Properties of the sampler, identical for both mixture summaries.",
    "\\item[c] Split-$\\hat{R}$ and effective sample size of the log posterior density (lp\\_\\_), which does not depend on component labels.",
    "\\item[d] Rank-normalized $\\hat{R}$ and bulk ESS of the component-specific parameters after Stephens relabeling within each fit, computed on 2,000 thinned draws; ESS is therefore bounded above by 2,000.",
    "\\item Single-component and frailty columns: split-$\\hat{R}$ and ESS over the intercept, coefficients, Weibull shape, global shrinkage scale and (frailty) $\\sigma_u$, taking the worst parameter in each fit.",
    "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
  )
}

#' Table 4: the largest between-cause differences in single-component coefficients, from paired stacked draws.
table_swings_tex <- function(sw) {
  labs <- c(design_labels(), clinical_labels())
  keep <- c("COHORT", "FACTOR_B", "ARM", "AMB_STAT3CAT", "GNRI")
  sw <- sw[match(keep, sw$variable), ]
  body <- vapply(seq_len(nrow(sw)), function(r) {
    z <- sw[r, ]
    sprintf(
      "%s & %s & %s & %s & %s \\\\",
      labs[z$variable], tex_cri(z$amp_median, z$amp_lo95, z$amp_hi95), tex_cri(z$death_median, z$death_lo95, z$death_hi95),
      tex_cri(z$diff_median, z$diff_lo95, z$diff_hi95, bold = excludes_zero(z$diff_lo95, z$diff_hi95)), tex_num(z$diff_sd, 2)
    )
  }, character(1))
  c(
    tex_header(),
    "\\begin{table}[htbp]", "\\centering", "\\footnotesize", "\\begin{threeparttable}",
    "\\caption{Largest between-cause differences in single-component coefficients, excluding age.}\\label{tab:swings}",
    "\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{@{}lcccc@{}}", "\\toprule",
    " & Amputation & Death & Amputation $-$ death & SD\\tnote{a} \\\\", "\\midrule", body,
    "\\bottomrule", "\\end{tabular}",
    "\\begin{tablenotes}[flushleft]\\scriptsize",
    "\\item[a] Posterior standard deviation of the difference.",
    "\\item Posterior median (95\\% equal-tailed credible interval) on the log-time scale. Cohort 2 = no adequate single-segment great saphenous vein (GSV) conduit; endovascular arm vs.\\ surgical bypass. Differences were computed draw by draw within each imputed dataset, pairing independent posterior draws of the two cause-specific fits, and then stacked across the 50 datasets; this retains the correlation the shared imputation induces between the two causes. Bold: interval excludes zero.",
    "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
  )
}
