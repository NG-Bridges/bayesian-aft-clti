#' Write a figure as PDF and PNG from one drawing function.
save_figure <- function(base, draw, width, height, res = 400, serif = TRUE) {
  fam <- if (serif) "serif" else "sans"
  pdf(paste0(base, ".pdf"), width = width, height = height, family = fam)
  draw()
  invisible(dev.off())
  pa <- list(filename = paste0(base, ".png"), width = width, height = height, units = "in", res = res, family = fam)
  if (capabilities("cairo")) pa$type <- "cairo"
  do.call(png, pa)
  draw()
  invisible(dev.off())
}

#' Stack the Phase 1 AFT_single draws per cause and check they reproduce the pooled quantiles.
stacked_single_draws <- function(ckpt_dir, co, n_imputations) {
  draws <- list()
  for (cs in names(cause_labels())) {
    files <- file.path(ckpt_dir, sprintf("inf_AFT_single_%s_%d.rds", cs, seq_len(n_imputations)))
    if (!all(file.exists(files))) stop("Missing Phase 1 checkpoints in ", ckpt_dir)
    d <- do.call(rbind, lapply(files, function(f) readRDS(f)$draws))
    z <- co[co$cause == cs & co$variable %in% colnames(d), ]
    q <- t(apply(d[, z$variable, drop = FALSE], 2, quantile, c(.5, .025, .975), names = FALSE))
    err <- max(abs(q - as.matrix(z[, c("post_median", "lo95", "hi95")])))
    if (err > 1e-6) stop("Checkpoint draws do not reproduce phase1_coefficients_pooled.csv (", cs, ", max diff ", signif(err, 3), ")")
    draws[[cs]] <- d
  }
  draws
}

#' Figure 1: AFT time ratios for stated clinical contrasts (A) and the eGFR time-ratio curve (B).
#' Penalized predictors are per SD, so log TR = beta * step / sd(x); eGFR_sq is the square of
#' standardized eGFR, so log TR(e vs r) = b1 (z_e - z_r) + b2 (z_e^2 - z_r^2).
figure_aft_forest <- function(co, mom, vi, draws, base) {
  cfg <- list(
    xlim = c(0.2, 2.6), ticks = c(0.25, 0.5, 1, 2),
    egfr_ref = 90, egfr_range = c(15, 120), # eGFR is winsorized to [5, 120]
    curve_draws = 20000, width = 7.5, height = 9.8, res = 400, serif = TRUE
  )
  red <- "#e34948"
  blue <- "#2a78d6"
  gray <- "#86857f"
  ink <- "#0b0b0b"
  ink2 <- "#52514e"
  band <- "#f1f0ec"
  curve_band <- "#d9d8d2"

  row_spec <- function(g, v, lab, con, step = 1) data.frame(group = g, var = v, lab = lab, con = con, step = step, stringsAsFactors = FALSE)
  # ESKD (kidney transplant, n = 21) is omitted from the figure; it is reported in the tables
  spec <- rbind(
    row_spec("Limb and anatomy", "Wound", "WIfI wound grade", "+1 grade"),
    row_spec("Limb and anatomy", "ISCHEMIA", "WIfI ischemia grade", "+1 grade"),
    row_spec("Limb and anatomy", "INFECTION", "WIfI foot infection grade", "+1 grade"),
    row_spec("Limb and anatomy", "FACTOR_B", "Infrapopliteal disease", "present"),
    row_spec("Limb and anatomy", "PREV_INFR_IND", "Prior infrainguinal revasc.", "yes vs. no"),
    row_spec("Cardiorenal and comorbidity", "eGFR", "eGFR", "see panel B"),
    row_spec("Cardiorenal and comorbidity", "prior_cardiac_dx", "Prior cardiac disease", "yes vs. no"),
    row_spec("Cardiorenal and comorbidity", "CHF_NYHA", "Congestive heart failure", "yes vs. no"),
    row_spec("Cardiorenal and comorbidity", "cerebrovascular_dx", "Cerebrovascular disease", "yes vs. no"),
    row_spec("Cardiorenal and comorbidity", "COPD_YN", "COPD", "yes vs. no"),
    row_spec("Cardiorenal and comorbidity", "DIAB_YN", "Diabetes mellitus", "yes vs. no"),
    row_spec("Cardiorenal and comorbidity", "SMOKE_HX", "Smoking history", "+1 level"),
    row_spec("Cardiorenal and comorbidity", "STATIN", "Statin use", "yes vs. no"),
    row_spec("Function, nutrition, demographics", "AGE_T", "Age", "+10 years"),
    row_spec("Function, nutrition, demographics", "GENDER", "Female sex", "vs. male"),
    row_spec("Function, nutrition, demographics", "AMB_STAT3CAT", "Ambulatory status", "1 level worse"),
    row_spec("Function, nutrition, demographics", "GNRI", "GNRI", "+10 points", 10),
    row_spec("Trial design", "COHORT", "No adequate GSV", "Cohort 2 vs. 1"),
    row_spec("Trial design", "ARM", "Endovascular therapy", "vs. bypass")
  )
  design_vars <- c("COHORT", "FACTOR_B", "ARM")
  cause <- cause_labels()
  msd <- setNames(mom$sd, mom$variable)
  mmu <- setNames(mom$mean, mom$variable)
  if (length(miss <- setdiff(setdiff(spec$var, design_vars), names(msd)))) stop("predictor_moments.csv lacks: ", paste(miss, collapse = ", "))
  if (length(miss <- setdiff(c(spec$var, "eGFR_sq"), co$variable))) stop("coefficients lack: ", paste(miss, collapse = ", "))

  tr_of <- function(v, cs, step) {
    k <- if (v %in% design_vars) 1 else step / msd[[v]]
    exp(quantile(draws[[cs]][, v] * k, c(.5, .025, .975), names = FALSE))
  }
  rows <- list()
  for (g in unique(spec$group)) {
    rows[[length(rows) + 1]] <- list(type = "header", lab = g)
    for (j in which(spec$group == g)) rows[[length(rows) + 1]] <- c(list(type = "var"), as.list(spec[j, ]))
  }
  n <- length(rows)
  yv <- rev(seq_len(n)) - 1
  est <- list()
  for (cs in names(cause)) {
    est[[cs]] <- lapply(rows, function(r) if (r$type == "var" && r$var != "eGFR") tr_of(r$var, cs, r$step) else NULL)
  }
  rank_of <- function(v, cs) {
    if (!v %in% vi$variable) {
      return(NA_integer_)
    }
    vi[vi$variable == v, if (cs == "amputation") "amp_rank" else "death_rank"]
  }

  eg <- seq(cfg$egfr_range[1], cfg$egfr_range[2], by = 1)
  zf <- function(x) (x - mmu[["eGFR"]]) / msd[["eGFR"]]
  dz <- zf(eg) - zf(cfg$egfr_ref)
  dz2 <- zf(eg)^2 - zf(cfg$egfr_ref)^2
  curve <- list()
  for (cs in names(cause)) {
    d <- draws[[cs]]
    idx <- thin_index(nrow(d), cfg$curve_draws)
    l <- outer(d[idx, "eGFR"], dz) + outer(d[idx, "eGFR_sq"], dz2)
    curve[[cs]] <- exp(t(apply(l, 2, quantile, c(.5, .025, .975), names = FALSE)))
  }

  lx <- function(p) 10^(par("usr")[1] + p * diff(par("usr")[1:2])) # axis fraction -> data (log x)
  y_at_line <- function(l) par("usr")[3] - l * par("csi") * diff(par("usr")[3:4]) / par("pin")[2]
  ylim_a <- c(-0.6, n - 0.4)
  zebra <- function() {
    k <- 0
    for (i in seq_len(n)) {
      if (rows[[i]]$type == "header") {
        k <- 0
        next
      }
      if (k %% 2 == 0) rect(par("usr")[1] - 1, yv[i] - 0.5, par("usr")[2] + 1, yv[i] + 0.5, col = band, border = NA)
      k <- k + 1
    }
  }

  draw <- function() {
    par(family = if (cfg$serif) "serif" else "sans", las = 1, xpd = FALSE)
    layout(rbind(c(1, 2, 3, 4, 5), c(6, 6, 6, 6, 6), c(7, 8, 9, 10, 11)),
      widths = c(2.45, 2.15, 0.42, 2.15, 0.42), heights = c(6.55, 0.42, 2.45)
    )
    m_a <- c(3.6, 0.3, 3.3, 0.3)

    par(mar = m_a)
    plot.new()
    plot.window(c(0, 1), ylim_a, yaxs = "i")
    zebra()
    for (i in seq_len(n)) {
      r <- rows[[i]]
      if (r$type == "header") {
        text(0, yv[i], r$lab, adj = c(0, 0.5), font = 2, cex = 0.86, col = ink, xpd = NA)
      } else {
        text(0.04, yv[i], r$lab, adj = c(0, 0.5), cex = 0.84, col = ink)
        text(1, yv[i], r$con, adj = c(1, 0.5), cex = 0.74, font = 3, col = ink2)
      }
    }
    mtext("A", side = 3, line = 1.6, adj = 0, font = 2, cex = 1.05)

    forest <- function(cs) {
      par(mar = m_a)
      plot(NA, log = "x", xlim = cfg$xlim, ylim = ylim_a, xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
      zebra()
      abline(v = 1, lty = 2, col = ink2, lwd = 0.8)
      for (i in seq_len(n)) {
        r <- rows[[i]]
        if (r$type != "var") next
        if (r$var == "eGFR") {
          lab <- "nonlinear: see panel B"
          w <- strwidth(lab, cex = 0.7, font = 3) / 2
          c0 <- 10^mean(par("usr")[1:2])
          rect(c0 / 10^w * 0.97, yv[i] - 0.32, c0 * 10^w * 1.03, yv[i] + 0.32, col = band, border = NA)
          text(c0, yv[i], lab, cex = 0.7, font = 3, col = ink2)
          next
        }
        e <- est[[cs]][[i]]
        m <- e[1]
        lo <- e[2]
        hi <- e[3]
        col <- if (hi < 1) red else if (lo > 1) blue else gray
        lo_c <- max(lo, cfg$xlim[1] * 1.01)
        hi_c <- min(hi, cfg$xlim[2] / 1.01)
        segments(lo_c, yv[i], hi_c, yv[i], col = col, lwd = 1.7, lend = 1)
        if (lo < cfg$xlim[1]) arrows(lo_c * 1.15, yv[i], lo_c, yv[i], length = 0.04, col = col, lwd = 1.4)
        if (hi > cfg$xlim[2]) arrows(hi_c / 1.15, yv[i], hi_c, yv[i], length = 0.04, col = col, lwd = 1.4)
        points(m, yv[i], pch = 21, cex = 1.05, lwd = 1.5, col = col, bg = if (col == gray) "white" else col)
      }
      axis(1, at = cfg$ticks, labels = as.character(cfg$ticks), cex.axis = 0.8, col = ink2, col.axis = ink, tcl = -0.3, mgp = c(2, 0.4, 0))
      mtext(cause[[cs]], side = 3, line = 1.6, font = 2, cex = 0.8, col = ink)
      mtext("Time ratio (95% CrI)", side = 3, line = 0.5, cex = 0.66, col = ink2)
      ya <- y_at_line(1.9)
      yt <- y_at_line(2.75)
      par(xpd = NA)
      arrows(lx(0.44), ya, lx(0.04), ya, length = 0.05, col = ink2, lwd = 0.8)
      arrows(lx(0.56), ya, lx(0.96), ya, length = 0.05, col = ink2, lwd = 0.8)
      text(lx(0.24), yt, "Earlier event", cex = 0.72, col = ink2)
      text(lx(0.76), yt, "Later event", cex = 0.72, col = ink2)
      par(xpd = FALSE)
    }
    rankcol <- function(cs) {
      par(mar = m_a)
      plot.new()
      plot.window(c(0, 1), ylim_a, yaxs = "i")
      zebra()
      for (i in seq_len(n)) {
        r <- rows[[i]]
        if (r$type != "var") next
        rk <- rank_of(r$var, cs)
        if (is.na(rk)) next
        text(0.5, yv[i], rk, cex = 0.76, font = if (rk <= 5) 2 else 1, col = if (rk <= 5) ink else ink2)
      }
      mtext("RSF\nrank", side = 3, line = 0.4, cex = 0.62, col = ink2)
    }
    forest("amputation")
    rankcol("amputation")
    forest("death")
    rankcol("death")

    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center",
      ncol = 3, bty = "n", cex = 0.76, x.intersp = 0.7,
      text.width = c(0.24, 0.24, 0.17),
      legend = c("Credibly shorter time to event", "Credibly longer time to event", "95% CrI includes 1"),
      col = c(red, blue, gray), pt.bg = c(red, blue, "white"), pch = 21, lty = 1, lwd = 1.7,
      pt.lwd = 1.5, pt.cex = 1.05, seg.len = 1.6, text.col = ink
    )

    m_b <- c(3.4, 0.3, 2.2, 0.3)
    par(mar = m_b)
    plot.new()
    plot.window(c(0, 1), c(0, 1))
    mtext("B", side = 3, line = 0.9, adj = 0, font = 2, cex = 1.05)
    text(0, 0.80, "eGFR", adj = c(0, 0.5), font = 2, cex = 0.86, col = ink)
    text(0.04, 0.62, sprintf("Time ratio vs. eGFR %g", cfg$egfr_ref), adj = c(0, 0.5), cex = 0.78, col = ink)
    text(0.04, 0.48, "Posterior median and", adj = c(0, 0.5), cex = 0.74, font = 3, col = ink2)
    text(0.04, 0.36, "95% credible band", adj = c(0, 0.5), cex = 0.74, font = 3, col = ink2)

    curve_panel <- function(cs) {
      par(mar = m_b)
      plot(NA, log = "y", xlim = cfg$egfr_range, ylim = cfg$xlim, xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
      cv <- curve[[cs]]
      polygon(c(eg, rev(eg)), c(cv[, 2], rev(cv[, 3])), col = curve_band, border = NA)
      abline(h = 1, lty = 2, col = ink2, lwd = 0.8)
      abline(v = cfg$egfr_ref, lty = 3, col = ink2, lwd = 0.8)
      lines(eg, cv[, 1], col = ink, lwd = 1.8)
      axis(1,
        at = c(15, 30, 60, 90, 120)[c(15, 30, 60, 90, 120) >= cfg$egfr_range[1]],
        cex.axis = 0.78, col = ink2, tcl = -0.3, mgp = c(2, 0.4, 0)
      )
      axis(2, at = cfg$ticks, labels = as.character(cfg$ticks), cex.axis = 0.75, col = ink2, tcl = -0.3, mgp = c(2, 0.4, 0), las = 1)
      mtext(expression("eGFR (mL/min/1.73 m"^2 * ")"), side = 1, line = 2.0, cex = 0.66, col = ink2)
      mtext(cause[[cs]], side = 3, line = 0.6, font = 2, cex = 0.72, col = ink)
      par(xpd = NA)
      text(cfg$egfr_ref, cfg$xlim[2], "ref.", adj = c(0.5, -0.3), cex = 0.62, col = ink2)
      par(xpd = FALSE)
    }
    curve_panel("amputation")
    par(mar = m_b)
    plot.new()
    curve_panel("death")
    par(mar = m_b)
    plot.new()
  }
  save_figure(base, draw, cfg$width, cfg$height, cfg$res, cfg$serif)
}

#' Figure 2: cross-validated AUC and IPA by horizon, single-component AFT vs RSF.
#' Error bars are SDs across imputed datasets (same folds in each), not CV standard errors.
figure_cv <- function(bay, rsf, base) {
  col <- c(AFT_single = "#B03A2E", AFT_frailty = "#7D3C98", AFT_mixture = "#1E8449", RSF = "#1F618D")
  pch <- c(AFT_single = 16, AFT_frailty = 17, AFT_mixture = 18, RSF = 15)
  nm <- c(AFT_single = "AFT", AFT_frailty = "Frailty AFT", AFT_mixture = "Gated mixture AFT", RSF = "Random survival forest")
  cause <- cause_labels()
  ebar_v <- function(x, y0, y1, col, lwd = 1.1) suppressWarnings(arrows(x, y0, x, y1, angle = 90, code = 3, length = 0.025, col = col, lwd = lwd))

  rsf$model <- "RSF"
  keep <- intersect(unique(bay$imputation), unique(rsf$imputation))
  cv <- rbind(bay[bay$imputation %in% keep, names(rsf)], rsf[rsf$imputation %in% keep, ])
  mods <- c(intersect(c("AFT_single", "AFT_frailty", "AFT_mixture"), unique(cv$model)), "RSF")
  mo <- sort(unique(cv$months))
  xs <- seq_along(mo)

  draw <- function() {
    par(mfrow = c(2, 2), mar = c(3.7, 4.2, 2.4, 0.8), mgp = c(2.4, 0.6, 0), cex.axis = 0.85, cex.lab = 0.95, las = 1)
    for (cs in names(cause)) {
      for (met in c("AUC", "IPA")) {
        sub <- cv[cv$cause == cs, ]
        a <- do.call(rbind, lapply(mods, function(m) {
          z <- sub[sub$model == m, ]
          d <- aggregate(z[[met]], by = list(months = z$months), FUN = function(v) c(m = mean(v, na.rm = TRUE), s = sd(v, na.rm = TRUE)))
          data.frame(model = m, months = d$months, mean = d$x[, "m"], sd = d$x[, "s"])
        }))
        yl <- range(c(a$mean - a$sd, a$mean + a$sd), na.rm = TRUE)
        yl <- yl + c(-1, 1) * 0.08 * diff(yl)
        ref <- if (met == "AUC") 0.5 else 0
        yl[1] <- min(yl[1], ref - 0.02 * diff(yl))
        plot(NA,
          xlim = c(0.8, length(mo) + 0.2), ylim = yl, xaxt = "n", xlab = "",
          ylab = if (met == "AUC") "Time-dependent AUC" else "Index of Prediction Accuracy",
          main = sprintf("%s: %s", cause[cs], met), cex.main = 0.95, font.main = 1, bty = "l"
        )
        axis(1, at = xs, labels = mo)
        mtext("Horizon (months)", side = 1, line = 2.1, cex = 0.75)
        abline(h = ref, col = "grey55", lty = 3)
        for (m in mods) {
          z <- a[a$model == m, ]
          z <- z[match(mo, z$months), ]
          if (all(is.na(z$mean))) next
          lines(xs, z$mean, col = col[m], lwd = 1.9)
          ebar_v(xs, z$mean - z$sd, z$mean + z$sd, col[m])
          points(xs, z$mean, col = col[m], pch = pch[m], cex = 0.95)
        }
        if (cs == names(cause)[1] && met == "AUC") {
          legend(0.75, yl[1] + 0.12 * diff(yl), nm[mods],
            yjust = 0, col = col[mods], pch = pch[mods],
            lwd = 1.9, bty = "n", cex = 0.72, seg.len = 1.2
          )
        }
      }
    }
  }
  save_figure(base, draw, width = 7.2, height = 5.8, res = 400, serif = TRUE)
}

#' Figure 3 and Figure S: calibration curves (natural-spline smooth of pseudo-values with a bootstrap band),
#' risk-decile points, and mirrored histograms of predicted risk under the axis.
figure_calibration <- function(summ, dec, crv, pbar, months, file, width, height, cfg) {
  mods <- cfg$calib_models
  pal <- c(AFT_single = "#2a78d6", RSF = "#eb6834", AFT_frailty = "#1baf7a")
  pch <- c(AFT_single = 16, RSF = 17, AFT_frailty = 15)
  lab <- c(AFT_single = "Single-component AFT", RSF = "Random survival forest", AFT_frailty = "Frailty AFT")
  ink <- "#1f1f1e"
  ink2 <- "#5f5e58"
  grid_col <- "#e4e3dd"
  alpha_hex <- function(col, a) rgb(t(col2rgb(col)) / 255, alpha = a)

  panel <- function(cs, mo) {
    h <- match(mo, cfg$horizons_months)
    cz <- crv[crv$cause == cs & crv$months == mo, ]
    dz <- dec[dec$cause == cs & dec$months == mo, ]
    sz <- summ[summ$cause == cs & summ$months == mo, ]
    pr <- unlist(lapply(mods, function(m) pbar[[m]][[cs]][, h]))
    lim <- min(1, 1.08 * max(quantile(pr, 0.99), dz$hi, cz$hi[is.finite(cz$hi)], na.rm = TRUE))
    strip <- 0.16 * lim
    plot(NULL, xlim = c(0, lim), ylim = c(-strip, lim), xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
    at <- pretty(c(0, lim), 4)
    at <- at[at <= lim]
    abline(v = at, h = at, col = grid_col, lwd = 0.8)
    axis(1, at = at, labels = sprintf("%.0f%%", 100 * at), col = ink2, col.axis = ink2, cex.axis = 0.8, lwd = 0.8)
    axis(2, at = at, labels = sprintf("%.0f%%", 100 * at), col = ink2, col.axis = ink2, cex.axis = 0.8, las = 1, lwd = 0.8)
    segments(0, 0, lim, lim, col = ink2, lty = 2, lwd = 1)
    br <- seq(0, lim, length.out = 41)
    mid <- -strip / 2
    hs <- lapply(mods, function(m) hist(pmin(pbar[[m]][[cs]][, h], lim), breaks = br, plot = FALSE)$counts)
    sc <- (strip / 2 - strip * 0.06) / max(unlist(hs))
    # First model's histogram points up, the second's down; there is no room for a third
    for (k in seq_along(mods)) {
      if (k > 2) next
      dirn <- if (k == 1) 1 else -1
      rect(br[-length(br)], mid, br[-1], mid + dirn * hs[[k]] * sc, col = alpha_hex(pal[mods[k]], 0.75), border = "white", lwd = 0.5)
    }
    segments(0, mid, lim, mid, col = ink2, lwd = 0.6)
    segments(0, 0, lim, 0, col = ink2, lwd = 0.8)
    for (m in mods) {
      z <- cz[cz$model == m, ]
      ok <- is.finite(z$lo) & is.finite(z$hi)
      polygon(c(z$predicted[ok], rev(z$predicted[ok])), c(z$lo[ok], rev(z$hi[ok])), col = alpha_hex(pal[m], 0.15), border = NA)
    }
    for (m in mods) {
      z <- cz[cz$model == m, ]
      lines(z$predicted, z$observed, col = pal[m], lwd = 2)
      q <- dz[dz$model == m, ]
      segments(q$predicted, pmax(q$lo, -strip * 0.02), q$predicted, pmin(q$hi, lim), col = pal[m], lwd = 1)
      points(q$predicted, q$observed, pch = pch[m], col = "white", cex = 1.25)
      points(q$predicted, q$observed, pch = pch[m], col = pal[m], cex = 0.95)
    }
    box(col = ink2, lwd = 0.8)
    title(main = sprintf("%s, %d months", tools::toTitleCase(cs), mo), cex.main = 0.95, col.main = ink, font.main = 2, line = 0.6)
    txt <- vapply(mods, function(m) {
      r <- sz[sz$model == m, ]
      sprintf("O/E %.2f  slope %.2f  ICI %.3f", r$OE, r$slope, r$ICI)
    }, character(1))
    legend("topleft", legend = txt, pch = pch[mods], col = pal[mods], text.col = ink, bty = "n", cex = 0.72, pt.cex = 0.95, inset = c(0.01, 0.01))
  }

  png(file, width = width, height = height, res = 150, bg = "white")
  nc <- length(months)
  layout(rbind(matrix(seq_len(2 * nc), 2, nc, byrow = TRUE), rep(2 * nc + 1, nc)), heights = c(1, 1, 0.14))
  par(mar = c(3.3, 3.9, 2.2, 0.8), mgp = c(2.2, 0.6, 0), tcl = -0.25, family = "sans")
  for (cs in names(cause_labels())) {
    for (mo in months) {
      panel(cs, mo)
      mtext("Predicted risk", side = 1, line = 2.1, cex = 0.7, col = ink2)
      mtext("Observed risk", side = 2, line = 2.7, cex = 0.7, col = ink2)
    }
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  labs <- c(lab[mods], "Ideal")
  legend("center",
    horiz = TRUE, bty = "n", cex = 0.85, text.col = ink, legend = labs,
    text.width = max(strwidth(labs, cex = 0.85)) * 1.05, x.intersp = 0.8,
    col = c(pal[mods], ink2), lty = c(rep(1, length(mods)), 2), lwd = c(rep(2, length(mods)), 1),
    pch = c(pch[mods], NA)
  )
  invisible(dev.off())
}
