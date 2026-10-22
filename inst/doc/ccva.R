## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4.4, dpi = 90, out.width = "100%")
set.seed(2026)

## -----------------------------------------------------------------------------
library(ccva)

## -----------------------------------------------------------------------------
d <- as_ccva(ccva_sim,
  exposure    = c("E1_risk", "E2_risk", "E3_risk"),
  sensitivity = grep("^s[0-9]+_", names(ccva_sim), value = TRUE),
  adaptive    = grep("^a[0-9]+_", names(ccva_sim), value = TRUE),
  outcomes    = c("iucn", "pop_trend", "range_change"),
  direction   = c(pop_trend = -1, range_change = -1),
  id = "species")
d$truth <- ccva_sim[c("E_star", "S_star", "A_star", "V_star")]
d

## -----------------------------------------------------------------------------
ccva_sdm_exposure(n_current = c(120, 80), n_future = c(60, 95),
                  n_overlap = c(50, 70), n_protected = c(30, 10))

## ----fig.height = 7.5---------------------------------------------------------
ccva_plot(d, "indicators")

## -----------------------------------------------------------------------------
chk <- ccva_check(d)
chk

## ----fig.height = 3-----------------------------------------------------------
ccva_plot(chk, "alpha")

## -----------------------------------------------------------------------------
ofori <- ccva_combine(d, "weighted", component_method = "ordinal")
ofori

## ----fig.height = 3.6---------------------------------------------------------
ccva_plot(ofori, "structure")

## ----fig.height = 3-----------------------------------------------------------
ccva_plot(ofori, "components")

## ----fig.height = 5.5---------------------------------------------------------
ccva_plot(ofori, "contributions", n = 8)

## -----------------------------------------------------------------------------
latent <- ccva_combine(d, "latent",
  shared_error = list(c("E1_risk", "E2_risk", "E3_risk")))
latent

## ----fig.height = 3.6---------------------------------------------------------
ccva_plot(latent, "structure")

## ----fig.height = 7.5---------------------------------------------------------
ccva_plot(latent, "effects")

## -----------------------------------------------------------------------------
cmp <- ccva_compare(ofori = ofori, latent = latent)
cmp

## ----fig.height = 4-----------------------------------------------------------
ccva_plot(cmp, "truth")

## ----fig.height = 5-----------------------------------------------------------
ccva_plot(cmp, "ranks")

## -----------------------------------------------------------------------------
ccva_default_constraints()

## -----------------------------------------------------------------------------
cons <- ccva_constraints(
  constrain_equal(c("a1_genetic", "a2_bodysize")),
  constrain_fixed("s7_repromode", 0),
  constrain_bounds("exposure", lower = 0, upper = 2))
cons
fit_c <- ccva_combine(d, "latent", constraints = cons,
  shared_error = list(c("E1_risk", "E2_risk", "E3_risk")))
fit_c$details$effects

## -----------------------------------------------------------------------------
small_n <- do.call(rbind, lapply(1:20, function(i) {
  s <- ccva_simulate(n = 50, seed = 100 + i)
  run <- function(cons) {
    f <- tryCatch(suppressWarnings(ccva_combine(s, "latent",
           constraints = cons,
           shared_error = list(c("E1_risk", "E2_risk", "E3_risk")))),
           error = function(e) NULL)
    if (is.null(f)) return(c(converged = 0, rho = NA))
    c(converged = f$details$converged,
      rho = cor(f$score, s$truth$V_star, method = "spearman"))
  }
  rbind(data.frame(model = "constrained", t(run(NULL))),
        data.frame(model = "free", t(run(ccva_unconstrained()))))
}))
aggregate(cbind(converged, rho) ~ model, small_n, mean, na.action = na.pass)

## ----fig.height = 3.6---------------------------------------------------------
mod <- ccva_combine(d, "shape",
                    constraints = ccva_constraints(constrain_moderator()))
ccva_plot(mod, "structure")

## ----fig.height = 3.2---------------------------------------------------------
shape <- ccva_combine(d, "shape")
ccva_plot(shape, "effects")

## -----------------------------------------------------------------------------
g <- ccva_simulate(n = 150, n_genomic = 300, n_causal = 10, seed = 3)
length(unlist(g$indicators))   # predictors
pen <- ccva_combine(g, "penalised")
pen

## ----fig.height = 5.5---------------------------------------------------------
ccva_plot(pen, "effects")

## -----------------------------------------------------------------------------
unc <- ccva_uncertainty(latent, "bootstrap", R = 40, seed = 1,
                        progress = FALSE)
unc

## ----fig.height = 6-----------------------------------------------------------
ccva_plot(unc, "intervals", n = 25)

## -----------------------------------------------------------------------------
ccva_uncertainty(ofori, "weights", R = 1000, seed = 1)

## -----------------------------------------------------------------------------
sapply(c("uniform", "quantile", "mixture", "roc"), function(m)
  table(ccva_classify(latent$score, m, data = d)))

## ----fig.height = 6-----------------------------------------------------------
ccva_plot(latent, "ranking", n = 8)

## -----------------------------------------------------------------------------
head(ccva_rank(latent, unc))

## -----------------------------------------------------------------------------
ccva_register("median", function(X, anchor, signs, ...) {
  z <- scale(X)
  ctr <- attr(z, "scaled:center"); scl <- attr(z, "scaled:scale")
  list(score = apply(z, 1, median),
       predict = function(M) apply(scale(M, ctr, scl), 1, median))
}, type = "component", description = "Median of standardised indicators")

ccva_sensitivity(d, method = "median")
ccva_unregister("median")

