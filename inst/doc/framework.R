## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4.4, dpi = 90, out.width = "100%")
set.seed(2026)

## -----------------------------------------------------------------------------
library(ccva)

## -----------------------------------------------------------------------------
ccva_recommend(ccva_simulate(n = 40, seed = 1))

## -----------------------------------------------------------------------------
regional <- ccva_simulate(n = 120, n_genomic = 250, seed = 2)
ccva_recommend(regional)

## -----------------------------------------------------------------------------
assays <- data.frame(species = rep(sprintf("sp_%02d", 1:25), each = 3),
                     ramp = rep(c(0.1, 0.25, 0.5), 25),
                     assay = rep(c("dynamic", "dynamic", "static"), 25))
true_ctmax <- rnorm(25, 44, 2)
assays$ctmax <- rep(true_ctmax, each = 3) + 3 * assays$ramp -
  1.2 * (assays$assay == "static") + rnorm(75, 0, 0.3)

adj <- ccva_adjust_method(assays, "ctmax", c("ramp", "assay"),
                          unit = "species")
adj
head(adj$by_unit)

## -----------------------------------------------------------------------------
d <- ccva_sim_data()
d$formative <- "adaptive"
ccva_check(d)

## -----------------------------------------------------------------------------
fw <- ccva_combine(d, "latent", formative_weights = "estimated",
                   shared_error = list(c("E1_risk", "E2_risk", "E3_risk")))
fw$details$indicator_weights

## -----------------------------------------------------------------------------
dm <- ccva_sim
dm$s1_thermal[sample(300, 40)] <- NA
dm$iucn[1:30] <- NA                       # e.g. Data Deficient species
xm <- as_ccva(dm, exposure = c("E1_risk", "E2_risk", "E3_risk"),
  sensitivity = grep("^s[0-9]+_", names(dm), value = TRUE),
  adaptive = grep("^a[0-9]+_", names(dm), value = TRUE),
  outcomes = c("iucn", "pop_trend", "range_change"),
  direction = c(pop_trend = -1, range_change = -1), id = "species",
  na_action = "keep")
fm <- ccva_combine(xm, "latent",
                   shared_error = list(c("E1_risk", "E2_risk", "E3_risk")))
cor(fm$score[1:30], ccva_sim$V_star[1:30], method = "spearman")

## -----------------------------------------------------------------------------
r <- ccva_combine(d, "rule")
r
head(r$details$drivers[order(-r$details$drivers$n_only), ])

## -----------------------------------------------------------------------------
s <- ccva_simulate(n = 400, link_AE = -0.6, seed = 4)
causal <- ccva_combine(s, "latent",
  shared_error = list(c("E1_risk", "E2_risk", "E3_risk")),
  constraints = ccva_constraints(
    constrain_link("adaptive", "exposure", "-"),
    constrain_link("adaptive", "sensitivity", "-")))
causal

## ----fig.height = 4-----------------------------------------------------------
ccva_plot(causal, "structure")

## -----------------------------------------------------------------------------
ccva_combine(s, "weighted", constraints = ccva_constraints(
  constrain_link("adaptive", "exposure")))

## ----fig.height = 3.2---------------------------------------------------------
nl <- ccva_simulate(n = 300, nonlinear = TRUE, seed = 7)
curved <- ccva_combine(nl, "shape", constraints = ccva_constraints(
  constrain_shape("exposure", "concave"),
  constrain_shape("sensitivity", "convex")))
ccva_plot(curved, "effects")

## -----------------------------------------------------------------------------
fit <- ccva_combine(d, "latent",
  shared_error = list(c("E1_risk", "E2_risk", "E3_risk")))
unc <- ccva_uncertainty(fit, "bootstrap", R = 40, seed = 1,
                        progress = FALSE)
risk <- ccva_calibrate(fit, outcome = "iucn", uncertainty = unc,
                       R = 200, seed = 1)
risk

## ----fig.height = 4-----------------------------------------------------------
ccva_plot(risk, "curve")

## ----fig.height = 5.5---------------------------------------------------------
ccva_plot(risk, "risk", n = 20)

## -----------------------------------------------------------------------------
cv <- ccva_cv(fit, k = 5, seed = 1, progress = FALSE)
cv
ccva_calibrate(fit, outcome = "iucn", cv = cv, R = 200, seed = 1)

## -----------------------------------------------------------------------------
occ <- unique(data.frame(site = rep(sprintf("site_%02d", 1:15), each = 30),
                         species = sample(names(fit$score), 450,
                                          replace = TRUE)))
sites <- ccva_aggregate(fit, occ, uncertainty = unc)
sites

## ----fig.height = 5-----------------------------------------------------------
ccva_plot(sites)

