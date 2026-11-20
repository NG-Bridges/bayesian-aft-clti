# Modeling decisions

## Competing risks, not illness-death

**Decision.** Amputation and death without amputation are modeled as cause-specific competing risks from randomization. Post-amputation death is not modeled.

**Why.** Only 93 patients died after an amputation, too few for a reliable parametric fit. A credible post-amputation model would also need updated clinical status after the amputation, and we only include covariates at baseline.

**Consequence.** Each cause is fit separately, with the competing event censored. Cumulative incidence is rebuilt from the pair of cause-specific survival curves.

## Weibull baseline hazard for both causes

**Decision.** Both cause-specific models use a Weibull baseline hazard.

**Why.**
- *Death without amputation.* The Weibull fits well. In the complete-case data, a generalized gamma did not improve on it (likelihood-ratio p = 0.54 and 0.17 in the two complete-case samples; Cox–Snell slopes ≈ 1). The estimated shape is close to 1 (1.10, 95% CI 1.02–1.19), so the hazard is close to constant. Adding a frailty leaves WAIC unchanged (ΔWAIC +1.0, SE 0.82).
- *Amputation, complete-case data.* The complete-case data favor the generalized gamma (ΔAIC 35–43; Cox–Snell slope 0.89 for the Weibull vs 1.00 for the generalized gamma).
- *Amputation, full imputed sample.* After imputing albumin for the full sample, the generalized gamma's shape parameter was no longer identified (Q = 1.50, 95% CI −6.7 to 9.7, fraction of missing information 0.99). Its location coefficients were nearly identical to the Weibull's. The data therefore could not support a claim about hazard shape, and the coefficients did not depend on it. The amputation shape is 0.58 (0.51–0.65), so the model already represents a falling hazard, with risk concentrated early after the procedure.

## Analysis cohort

The raw export has 1,832 records for 1,830 patients. Cleaning removes:
- 14 near-empty records (4 or more comorbidity fields missing)
- the 4 records of 2 patients who received both treatments (the export holds one record per index procedure,
  so these are the only patients with two records; both records carry the same randomized arm)

That leaves 1,814 patients. The analysis cohort is the 1,672 of them with every core predictor and a complete competing-risks outcome.

The 142 patients excluded at that last step are mostly missing a WIfI grade (ischemia 105, wound 58, foot infection 41) or eGFR (23).

## Multiple imputation of albumin, weight and height

**Decision.**
- Impute albumin, weight and height with predictive mean matching (`mice`, m = 50, 5 iterations), then compute GNRI in each completed dataset.
- Include in the imputation model the Nelson–Aalen cumulative hazards and event indicators for reintervention, amputation and death (White and Royston 2009).
- m = 50 follows the m ≥ 100 × FMI rule. GNRI's fraction of missing information is 0.50 for amputation and 0.37 for death. Every other coefficient is below 0.10.

**Consequence.** The imputation is run once, on the full sample, not inside each cross-validation training fold. We checked whether this leaks outcome information into the held-out folds. GNRI predicts no better in patients with imputed albumin than in those with measured albumin (the difference-in-differences in its AUC contribution is negative). A Cox interaction between GNRI and an imputed-albumin indicator is null for both causes (HR 0.90, p = 0.57; HR 1.15, p = 0.24).

## Albumin missingness is not a predictor

Patients without an albumin measurement have a lower hazard of death (HR 0.73, 0.60–0.89), presumably because labs were ordered for sicker-looking patients. Missingness is therefore never used as a model feature. If a reviewer raises missing-not-at-random concerns, a delta-adjusted sensitivity analysis is the planned response.

## Quadratic eGFR

eGFR enters with linear and quadratic terms. The quadratic term is the square of *standardized* eGFR, which keeps both terms on a similar scale. Curvature was significant for both causes (p = 0.008 for amputation, 1.1 × 10⁻⁶ for death), and a cubic spline added nothing beyond the quadratic. The final functional-form tests use a truncated power basis, so the linear, quadratic and spline models are properly nested.

## Standardization

- **Inference fits** use one fixed set of means and SDs, computed from the 50 completed datasets stacked together. Every fit is then on the same scale, and draws can be pooled directly.
- **Cross-validation fits** standardize on the training fold only.
- **Design variables** are coded 0/1 and never standardized.

## Mixture labels

**Decision.** The two mixture components are left unconstrained. Labels are resolved after fitting in two steps:
1. Stephens' algorithm, within each fit.
2. A cross-fit anchor: components are ordered by their allocation-weighted mean eGFR, so class 1 is the low-eGFR class in every dataset.

## Zero follow-up times

Twelve patients have a recorded event or censoring time of 0. These times are floored at half a day (0.5 / 365.25 years), because a log-time model is undefined at t = 0.
