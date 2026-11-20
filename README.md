# Accelerated Failure Time Models for Identifying Outcome-Specific Risk Factors and Assessing Latent Heterogeneity in Chronic Limb-Threatening Ischemia

Analysis code for Bridges N, McGinigle K, Hosseinian SM, "Accelerated Failure Time Models for Identifying Outcome-Specific Risk Factors and Assessing Latent Heterogeneity in Chronic Limb-Threatening Ischemia." The scripts reproduce the figures and tables in the article and its supplement from the BEST-CLI trial data.

The repository contains no patient-level data. The example run simulates a dataset of 400 patients with the same columns and coding as the trial export, so the code can be run end to end without the trial data.

## Setup

You need R 4.6 and rstan 2.26 or later. On Windows, also install the Rtools release that matches your R version. Install the pinned package versions with:

```sh
Rscript -e "renv::restore()"
```

## Running the example

```sh
CLTI_EXAMPLE=1 Rscript scripts/00_run_all.R
```

This runs the whole pipeline on the simulated data with shortened settings and writes the figures and tables to `output/example/`. A full example run takes about 80 minutes on an Apple M1 laptop. The example only checks that the code runs; its estimates are not interpretable.

## Running on the trial data

Place the BEST-CLI SAS datasets (the `.sas7bdat` files) in `data/raw/sas/` and run:

```sh
Rscript scripts/00_run_all.R
```

The full pipeline will likely take days to weeks, depending on the hardware. If a run is interrupted, rerunning it resumes from the last saved model fit.
