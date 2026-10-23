# ccva 0.1.1

## New features

* `ccva_combine(method = "rule")` implements the rule-based framework of
  Foden et al. (2013): indicators flag units as high exposure, high
  sensitivity or low adaptive capacity (coded 0/1 indicators at their risk
  end, continuous ones by a `threshold` quantile), a component is flagged if
  any indicator is, and units get Foden's categories (highly vulnerable,
  high latent risk, potential adapter, potential persister). Missing
  indicators are treated as `unknown = "optimistic"` or `"pessimistic"`.
  `details$drivers` and `ccva_plot(fit, "effects")` show how many units each
  indicator flags on its own.
* Missing values: `as_ccva(na_action = "keep")` keeps incomplete units. The
  `"mean"` and `"ordinal"` scorers average the observed indicators, and the
  latent model uses full-information maximum likelihood, so units with a
  missing outcome (e.g. IUCN Data Deficient species) are scored from their
  indicators. Diagnostics, validation and calibration use the available
  data; `summary()` and `print()` report missingness.
* Formative components with estimated weights:
  `ccva_combine(method = "latent", formative_weights = "estimated")` enters
  every indicator of a formative component as a cause of latent
  vulnerability (a MIMIC model) with its own sign-constrained weight,
  reported in `details$indicator_weights` and plotted by
  `ccva_plot(fit, "effects")`.
* Ordered categorical outcomes: `ccva_combine(method = "latent", ordered =)`
  fits outcomes such as IUCN categories as ordinal (WLSMV).
* `ccva_cv()` computes cross-validated scores: each fold is scored by a
  model refitted without it, with the outcomes of the held-out units hidden
  (as for an unassessed species). With `cluster`, whole groups (families,
  basins) are kept in one fold. `ccva_calibrate(cv = )` uses these scores
  for out-of-sample calibrated risk, and returns probability `draws` for
  intervals on sums such as the expected number of threatened species.
* Cluster bootstrap: `ccva_uncertainty(cluster = )` resamples whole groups
  of related or co-occurring units instead of single units, so intervals
  respect dependence within groups.
* Bootstrap and cross-validation refits of the latent model skip fit
  measures, which makes them several times faster with missing data.

## Bug fixes

* Shape-constrained fits now compute vulnerability scores with `predict()`
  on the fitted model rather than from the fitted values stored by `scam`,
  so that scores and `predict()` on the same data agree on every platform
  (fixes a test failure on r-oldrel-macos-arm64).

# ccva 0.1.0

First release.

* Data contract `as_ccva()`: any dataset, columns tagged to exposure,
  sensitivity and adaptive capacity, with directions and outcome proxies.
* Simulator `ccva_simulate()` with a known latent structure (linear or
  nonlinear truth, optional genomic block) and bundled dataset `ccva_sim`.
* Component scoring: Ofori-style ordinal, continuous mean, PCA and
  sign-constrained penalised regression.
* Combination models: weighted index (Ofori et al. 2017, corrected
  normalisation), shape-constrained additive model (scam), constrained
  latent SEM (lavaan) with shared-error method factors, and penalised
  two-stage model (glmnet).
* Constraint system on predictors and components: sign, bounds, equality,
  fixed values and adaptive capacity as moderator.
* Diagnostics, validation, method comparison, bootstrap and weight
  uncertainty (with factor-score error for latent models), calibrated
  classification, plots for every step and an HTML report.
* Method registry `ccva_register()` for contributed methods.
* Framework tools: directed links between components (`constrain_link()`)
  with direct, indirect and total effects, or two-stage resolution;
  formative components (`as_ccva(formative = )`) and mean inter-item
  correlation in `ccva_check()`; concave / convex shapes
  (`constrain_shape()`); calibrated risk predictions (`ccva_calibrate()`);
  method recommendations by data availability (`ccva_recommend()`); assay
  methodology adjustment (`ccva_adjust_method()`); community / site
  aggregation (`ccva_aggregate()`); workflow vignette "framework".
