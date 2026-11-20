#' The model specification that keys the saved fits.
#' Covers everything stored inside a checkpoint; grid_n is applied at scoring and is deliberately absent.
spec_string <- function(ms) {
  cfg <- ms$cfg
  paste(c(
    sort(ms$model_pen_vars), "|", sort(ms$design_vars), "|", sort(ms$gate_vars), "|",
    cfg$stan_warmup, cfg$stan_sample, cfg$stan_chains, cfg$K_mix,
    cfg$gh_nodes, cfg$design_scale, cfg$slab_scale, cfg$slab_df,
    cfg$p0_relevant, cfg$frailty_scale, cfg$intercept_scale, cfg$shape_scale,
    cfg$store_draws, cfg$pred_draws, cfg$relabel_draws, cfg$anchor_var
  ), collapse = "_")
}

spec_hash <- function(spec) {
  z <- utf8ToInt(spec)
  sprintf("%08x", as.integer(sum(z * seq_along(z)) %% 2147483647L))
}

checkpoint_dir <- function(ms, root = results_path("checkpoints")) {
  spec <- spec_string(ms)
  dir <- file.path(root, paste0("spec_", spec_hash(spec)))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  writeLines(spec, file.path(dir, "SPEC.txt"))
  # Lets later scripts find the fits without rebuilding the specification from patient data
  writeLines(basename(dir), file.path(root, "CURRENT"))
  dir
}

current_checkpoint_dir <- function(root = results_path("checkpoints")) {
  f <- file.path(root, "CURRENT")
  if (!file.exists(f)) stop("Missing ", f, "; run scripts/02_fit_models.R first.")
  file.path(root, readLines(f, n = 1))
}

checkpoint_path <- function(dir, ...) file.path(dir, paste0(paste(..., sep = "_"), ".rds"))

#' Return the saved result at `path`, or run `fn()` and save it.
#' With CLTI_NO_SAMPLING=1 a missing checkpoint is an error, so a check run never starts MCMC.
checkpointed <- function(path, fn) {
  if (file.exists(path)) {
    r <- try(readRDS(path), silent = TRUE)
    if (!inherits(r, "try-error")) {
      return(r)
    }
  }
  if (identical(Sys.getenv("CLTI_NO_SAMPLING"), "1")) {
    stop(structure(
      class = c("missing_checkpoint", "error", "condition"),
      list(message = paste("CLTI_NO_SAMPLING=1 and no checkpoint at", path), call = NULL)
    ))
  }
  r <- fn()
  saveRDS(r, path)
  r
}
