#' Cross-validated vulnerability scores
#'
#' The latent, shape-constrained and penalised models use the outcomes to
#' estimate vulnerability, so agreement between their scores and those same
#' outcomes is optimistic. `ccva_cv()` gives honest, out-of-sample scores:
#' the units are split into `k` folds; for each fold the whole pipeline is
#' refitted on the other folds, and the held-out units are scored with
#' their outcomes hidden -- exactly the position of a species that has never
#' been assessed, such as an IUCN Data Deficient species. The
#' cross-validated scores can then be compared with the hidden outcomes, or
#' passed to [ccva_calibrate()] for calibrated, out-of-sample risk.
#'
#' Related units (species of one family, or of one lake) are not
#' independent, and a random split lets a model learn a group from its
#' relatives. With `cluster`, whole groups are assigned to folds, so every
#' unit is scored by a model that has seen none of its group.
#'
#' The weighted index and the rule-based method do not use the outcomes,
#' so for them cross-validated and in-sample scores coincide up to the
#' data-dependent scoring thresholds; they are included so that all methods
#' can be compared in the same way.
#'
#' @param fit A `ccva_fit`.
#' @param k Number of folds.
#' @param hide Outcomes hidden for the held-out units; default all
#'   outcomes. Hiding only some (e.g. Red List status) mirrors units for
#'   which the other outcomes are known.
#' @param cluster Optional grouping of the units (a column name or one value
#'   per unit); groups are kept together when assigning folds.
#' @param seed Optional random seed.
#' @param progress Show progress.
#' @return An object of class `ccva_cv`: `score` (cross-validated scores on
#'   the scale of `fit$score`), `raw` (on the model's own scale), `fold`,
#'   and `validation` (Spearman correlation of in-sample and
#'   cross-validated scores with every hidden outcome).
#' @seealso [ccva_calibrate()], [ccva_validate()]
#' @export
#' @examples
#' \donttest{
#' d <- ccva_sim_data()
#' fit <- ccva_combine(d, "latent",
#'   shared_error = list(c("E1_risk", "E2_risk", "E3_risk")))
#' cv <- ccva_cv(fit, k = 5, seed = 1)
#' cv
#' # out-of-sample calibrated risk
#' ccva_calibrate(fit, outcome = "iucn", cv = cv, R = 100, seed = 1)
#' }
ccva_cv <- function(fit, k = 5, hide = NULL, cluster = NULL, seed = NULL,
                    progress = interactive()) {
  if (!inherits(fit, "ccva_fit")) stop("`fit` must be a ccva_fit.",
                                       call. = FALSE)
  x <- fit$data
  n <- length(x$id)
  if (!is.numeric(k) || length(k) != 1 || k < 2 || k > n)
    stop("`k` must be a number between 2 and the number of units.",
         call. = FALSE)
  hide <- hide %||% x$outcomes
  bad <- setdiff(hide, x$outcomes)
  if (length(bad))
    stop("`hide` must name outcomes; not outcomes: ",
         paste(bad, collapse = ", "), call. = FALSE)
  if (!is.null(fit$args$ordered) && length(intersect(hide, fit$args$ordered)))
    stop("Hidden outcomes cannot be ordered: scores with missing ordered ",
         "outcomes are not available. Refit with numeric outcomes.",
         call. = FALSE)
  if (!is.null(seed)) set.seed(seed)

  groups <- .cluster_groups(x, cluster)
  fold <- if (is.null(groups)) {
    sample(rep(seq_len(k), length.out = n))
  } else {
    ug <- unique(groups)
    if (length(ug) < k)
      stop("Fewer groups than folds; reduce `k`.", call. = FALSE)
    gf <- stats::setNames(sample(rep(seq_len(k), length.out = length(ug))),
                          ug)
    unname(gf[groups])
  }

  op <- options(ccva.light_refit = TRUE, ccva.force_fiml = TRUE)
  on.exit(options(op), add = TRUE)
  raw <- rep(NA_real_, n)
  for (f in seq_len(k)) {
    if (progress) message("fold ", f, " of ", k)
    train <- .subset_ccva(x, which(fold != f))
    fb <- suppressMessages(suppressWarnings(
      do.call(ccva_combine, c(list(x = train), fit$args))))
    nd <- x$data[fold == f, , drop = FALSE]
    nd[hide] <- NA
    raw[fold == f] <- fb$predict_raw(nd)
  }
  score <- stats::setNames(.rescale01(raw, fit$raw_range[1], fit$raw_range[2]),
                           x$id)

  O <- .oriented_outcomes(x)
  validation <- do.call(rbind, lapply(hide, function(o) {
    y <- O[, o]
    data.frame(outcome = o,
               in_sample = stats::cor(fit$score, y, method = "spearman",
                                      use = "complete.obs"),
               cross_validated = stats::cor(score, y, method = "spearman",
                                            use = "complete.obs"),
               n = sum(!is.na(y)), row.names = NULL)
  }))
  structure(list(score = score, raw = stats::setNames(raw, x$id),
                 fold = fold, k = k, hide = hide, validation = validation,
                 clustered = !is.null(groups), method = fit$method,
                 uses_outcomes = fit$method %in% c("latent", "shape",
                                                   "penalised")),
            class = "ccva_cv")
}

#' @export
print.ccva_cv <- function(x, ...) {
  cat(sprintf("<ccva_cv> %d-fold cross-validation%s of the %s model\n",
              x$k, if (x$clustered) " (whole groups per fold)" else "",
              x$method))
  cat("  hidden for held-out units:", paste(x$hide, collapse = ", "), "\n")
  if (!x$uses_outcomes)
    cat("  this method does not use the outcomes, so in-sample and",
        "cross-validated scores agree closely\n")
  cat("  Spearman correlation with each hidden outcome:\n")
  v <- x$validation
  v$in_sample <- round(v$in_sample, 3)
  v$cross_validated <- round(v$cross_validated, 3)
  print(v, row.names = FALSE)
  invisible(x)
}
