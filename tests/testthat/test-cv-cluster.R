test_that("ccva_cv gives out-of-sample scores for every unit", {
  d <- ccva_sim_data()
  se <- list(c("E1_risk", "E2_risk", "E3_risk"))
  fit <- ccva_combine(d, "latent", shared_error = se)
  cv <- ccva_cv(fit, k = 5, seed = 1, progress = FALSE)
  expect_s3_class(cv, "ccva_cv")
  expect_length(cv$score, length(d$id))
  expect_false(anyNA(cv$score))
  expect_setequal(unique(cv$fold), 1:5)
  v <- cv$validation
  expect_equal(nrow(v), 3)
  # hiding the outcomes costs agreement with them
  expect_true(all(v$cross_validated < v$in_sample))
  expect_gt(cor(cv$score, d$truth$V_star, method = "spearman"), 0.6)
  expect_output(print(cv), "cross-validation")
  expect_s3_class(ccva_plot(cv), "ggplot")
  expect_s3_class(ccva_plot(cv, "scores"), "ggplot")
  # only one outcome hidden
  cv1 <- ccva_cv(fit, k = 5, hide = "iucn", seed = 1, progress = FALSE)
  expect_equal(cv1$validation$outcome, "iucn")
  expect_error(ccva_cv(fit, hide = "E1_risk"), "outcomes")
  expect_error(ccva_cv(fit, k = 1), "between 2")
})

test_that("ccva_cv keeps groups together and works for other methods", {
  d <- ccva_sim_data()
  grp <- rep(seq_len(60), each = 5)
  w <- ccva_combine(d, "weighted")
  cvw <- ccva_cv(w, k = 4, cluster = grp, seed = 2, progress = FALSE)
  expect_true(cvw$clustered)
  expect_true(all(tapply(cvw$fold, grp, function(f) length(unique(f))) == 1))
  expect_gt(cor(cvw$score, w$score), 0.95)
  expect_error(ccva_cv(w, k = 4, cluster = rep(1:2, 150)), "at least 5")
})

test_that("calibration can use cross-validated scores", {
  d <- ccva_sim_data()
  fit <- ccva_combine(d, "latent",
                      shared_error = list(c("E1_risk", "E2_risk", "E3_risk")))
  cv <- ccva_cv(fit, k = 5, seed = 1, progress = FALSE)
  cal_in <- ccva_calibrate(fit, outcome = "iucn", R = 30, seed = 1)
  cal_cv <- ccva_calibrate(fit, outcome = "iucn", cv = cv, R = 30, seed = 1)
  expect_true(cal_in$in_sample)
  expect_false(cal_cv$in_sample)
  expect_true(cal_cv$cross_validated)
  expect_lt(cal_cv$metrics[["auc"]], cal_in$metrics[["auc"]])
  expect_equal(ncol(cal_cv$draws), length(d$id))
  expect_output(print(cal_cv), "cross-validated")
  expect_error(ccva_calibrate(fit, cv = list()), "ccva_cv")
})

test_that("cluster bootstrap resamples whole groups", {
  d <- ccva_sim_data()
  fit <- ccva_combine(d, "weighted")
  d$data$grp <- rep(seq_len(50), each = 6)
  fit$data <- d
  u <- ccva_uncertainty(fit, "bootstrap", R = 20, seed = 1, cluster = "grp",
                        progress = FALSE)
  expect_equal(u$cluster$n_groups, 50)
  expect_equal(u$cluster$name, "grp")
  expect_output(print(u), "cluster bootstrap over 50 groups")
  set.seed(3)
  rows <- ccva:::.cluster_rows(rep(1:10, each = 3))
  expect_equal(length(rows) %% 3, 0)
  expect_error(ccva_uncertainty(fit, "bootstrap", cluster = "nope"),
               "not found")
  expect_error(ccva_uncertainty(fit, "weights", cluster = "grp"),
               "bootstrap")
})
