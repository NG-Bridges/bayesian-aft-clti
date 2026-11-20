#' Competing-risks random survival forest on raw predictors (trees are scale-invariant; no eGFR_sq).
#' @param vars Predictors; the leakage check drops GNRI.
rsf_fit <- function(tr, importance, ms, vars = c(ms$pen_vars, ms$design_vars)) {
  cfg <- ms$cfg
  form <- as.formula(paste0("Surv(cr_time, cr_status) ~ ", paste(vars, collapse = " + ")))
  randomForestSRC::rfsrc(form,
    data = tr[, c("cr_time", "cr_status", vars)],
    ntree = cfg$rsf_ntree, nodesize = cfg$rsf_nodesize, nsplit = cfg$rsf_nsplit,
    splitrule = "logrankCR", importance = importance, seed = -abs(cfg$seed)
  )
}

#' RSF cumulative incidence at the horizons: the last forest time at or before each one.
rsf_cif <- function(fit, te, ms, vars = c(ms$pen_vars, ms$design_vars)) {
  pr <- predict(fit, newdata = te[, vars, drop = FALSE])
  if (is.null(pr$cif)) stop("rfsrc returned no CIF; cr_status must be 0/1/2 with both events in the fold.")
  ix <- sapply(ms$cfg$horizons, function(h) {
    j <- which(pr$time.interest <= h)
    if (length(j)) max(j) else 1L
  })
  nr <- dim(pr$cif)[1]
  list(amputation = matrix(pr$cif[, ix, 1], nrow = nr), death = matrix(pr$cif[, ix, 2], nrow = nr))
}

#' Permutation importance for both causes (rfsrc orders the columns by event code).
rsf_importance <- function(d, ms) {
  im <- rsf_fit(d, "permute", ms)$importance
  if (!is.matrix(im) || ncol(im) < 2) stop("rfsrc importance is not a competing-risks matrix.")
  list(
    tab = data.frame(variable = rownames(im), amputation = as.numeric(im[, 1]), death = as.numeric(im[, 2]), row.names = NULL),
    colnames = colnames(im)
  )
}

#' Out-of-fold RSF CIFs for one completed dataset.
rsf_oof <- function(d, fold_id, ms, vars = c(ms$pen_vars, ms$design_vars)) {
  h <- length(ms$cfg$horizons)
  a <- matrix(NA_real_, nrow(d), h)
  dm <- a
  for (k in seq_len(ms$cfg$n_folds)) {
    pr <- rsf_cif(rsf_fit(d[fold_id != k, ], "none", ms, vars), d[fold_id == k, ], ms, vars)
    a[fold_id == k, ] <- pr$amputation
    dm[fold_id == k, ] <- pr$death
  }
  list(amputation = a, death = dm)
}
