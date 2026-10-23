sim_with_missing <- function(n = 250, seed = 21, p_miss = 0.1) {
  s <- ccva_simulate(n = n, seed = seed)
  dat <- s$data
  set.seed(seed)
  ind <- unlist(s$indicators)
  for (v in ind) dat[[v]][stats::runif(n) < p_miss] <- NA
  dat$iucn[1:40] <- NA                      # e.g. Data Deficient
  list(dat = dat, sim = s)
}

build <- function(dat, s, ...) {
  as_ccva(dat, exposure = s$indicators$exposure,
          sensitivity = s$indicators$sensitivity,
          adaptive = s$indicators$adaptive,
          outcomes = c("iucn", "pop_trend", "range_change"),
          direction = c(pop_trend = -1, range_change = -1),
          id = "species", ...)
}

test_that("na_action = 'keep' retains units and records missingness", {
  m <- sim_with_missing()
  expect_message(x_omit <- build(m$dat, m$sim), "Dropping")
  x <- build(m$dat, m$sim, na_action = "keep")
  expect_gt(length(x$id), length(x_omit$id))
  expect_true(x$has_missing)
  expect_output(print(x), "Missing values")
  expect_true("n_missing" %in% names(summary(x)))
})

test_that("mean scorer and weighted index handle missing indicators", {
  m <- sim_with_missing()
  x <- build(m$dat, m$sim, na_action = "keep")
  e <- ccva_exposure(x)
  expect_false(anyNA(e$score))
  w <- ccva_combine(x, "weighted")
  expect_false(anyNA(w$score))
  expect_error(ccva_combine(x, "penalised"), "missing values")
  chk <- ccva_check(x)
  expect_true(all(is.finite(chk$alpha$alpha)))
})

test_that("latent model uses FIML and scores units with missing outcomes", {
  m <- sim_with_missing()
  x <- build(m$dat, m$sim, na_action = "keep")
  se <- list(c("E1_risk", "E2_risk", "E3_risk"))
  f <- ccva_combine(x, "latent", shared_error = se)
  expect_true(f$details$converged)
  expect_true(f$details$missing %in% c("ml", "ml.x"))
  expect_false(anyNA(f$score))
  expect_length(f$score, length(x$id))
  truth <- m$sim$truth$V_star[match(x$id, m$sim$data$species)]
  expect_gt(cor(f$score, truth, method = "spearman"), 0.8)
  # units whose IUCN status is missing are still scored sensibly
  miss <- is.na(x$data$iucn)
  expect_gt(cor(f$score[miss], truth[miss], method = "spearman"), 0.7)
  v <- ccva_validate(f)
  expect_true(all(is.finite(v$spearman)))
  r <- ccva_calibrate(f, outcome = "iucn", R = 20, seed = 1)
  expect_equal(nrow(r$units), length(x$id))
  expect_true(is.finite(r$metrics[["auc"]]))
})

test_that("formative components with estimated weights respect signs", {
  d <- ccva_sim_data()
  x <- d
  x$formative <- c("sensitivity", "adaptive")
  f <- ccva_combine(x, "latent", formative_weights = "estimated",
                    shared_error = list(c("E1_risk", "E2_risk", "E3_risk")))
  expect_true(f$details$converged)
  w <- f$details$indicator_weights
  expect_setequal(unique(w$component), c("sensitivity", "adaptive"))
  expect_true(all(w$estimate[w$component == "sensitivity"] >= -1e-6))
  expect_true(all(w$estimate[w$component == "adaptive"] <= 1e-6))
  e <- f$details$effects
  expect_equal(nrow(e), 3)
  expect_gt(e$std[e$component == "sensitivity"], 0)
  expect_lt(e$std[e$component == "adaptive"], 0)
  expect_gt(cor(f$score, d$truth$V_star, method = "spearman"), 0.8)
  expect_s3_class(ccva_plot(f, "effects"), "ggplot")
  expect_equal(unname(predict(f, d)), unname(f$score), tolerance = 1e-6)
  expect_error(ccva_combine(x, "latent", formative_weights = "estimated",
                            constraints = ccva_constraints(
                              constrain_link("adaptive", "sensitivity"))),
               "formative")
})

test_that("ordered outcomes are fitted with WLSMV on complete data", {
  d <- ccva_sim_data()
  f <- ccva_combine(d, "latent", ordered = c("iucn", "pop_trend"),
                    shared_error = list(c("E1_risk", "E2_risk", "E3_risk")))
  expect_true(f$details$converged)
  expect_equal(lavaan::lavInspect(f$model, "options")$estimator, "DWLS")
  expect_gt(cor(f$score, d$truth$V_star, method = "spearman"), 0.8)
  expect_error(ccva_combine(d, "latent", ordered = "E1_risk"), "outcome")
  m <- sim_with_missing()
  x <- build(m$dat, m$sim, na_action = "keep")
  expect_error(ccva_combine(x, "latent", ordered = "iucn"), "complete data")
})

test_that("rule method reproduces Foden-style classification", {
  dat <- data.frame(
    id = paste0("u", 1:8),
    e1 = c(1, 1, 1, 0, 1, 0, 0, NA),
    e2 = c(0, 0, 0, 0, 0, 0, 0, 0),
    s1 = c(1, 1, 0, 1, 1, 0, 1, 1),
    s2 = c(0, NA, 0, 0, 0, 0, 0, 0),
    a1 = c(1, 1, 1, 1, 0, 0, 1, 1),   # 1 = poor disperser (low capacity)
    y = 1:8)
  x <- as_ccva(dat, exposure = c("e1", "e2"), sensitivity = c("s1", "s2"),
               adaptive = "a1", outcomes = "y", id = "id",
               direction = c(a1 = -1), na_action = "keep")
  f <- ccva_combine(x, "rule")
  cat <- as.character(f$details$category)
  expect_equal(cat[1:4], c("highly vulnerable", "highly vulnerable",
                           "potential persister", "high latent risk"))
  expect_equal(cat[5], "potential adapter")
  expect_equal(cat[6], "lower risk")
  # unknown exposure: optimistic = not flagged, pessimistic = flagged
  expect_equal(cat[8], "high latent risk")
  fp <- ccva_combine(x, "rule", unknown = "pessimistic")
  expect_equal(as.character(fp$details$category)[8], "highly vulnerable")
  expect_equal(as.character(f$class)[1], "High")
  expect_equal(unname(f$score[1]), 1)
  expect_equal(unname(predict(f, x)), unname(f$score))
  dr <- f$details$drivers
  expect_equal(dr$n_only[dr$indicator == "s1"], 6)
  expect_output(print(f), "highly vulnerable")
  expect_s3_class(ccva_plot(f, "effects"), "ggplot")
  expect_s3_class(ccva_plot(f, "structure"), "ggplot")
})

test_that("rule method uses quantile thresholds for continuous indicators", {
  d <- ccva_sim_data()
  f <- ccva_combine(d, "rule", threshold = 0.25)
  dr <- f$details$drivers
  expect_true(all(dr$type == "quantile"))
  expect_true(all(abs(dr$n_flagged / length(d$id) - 0.25) < 0.02))
  expect_gt(cor(f$score, d$truth$V_star, method = "spearman"), 0.3)
  u <- ccva_uncertainty(f, R = 15, seed = 1, progress = FALSE)
  expect_s3_class(u, "ccva_uncertainty")
})
