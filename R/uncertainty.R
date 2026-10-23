#' Quantify uncertainty in vulnerability scores and rankings
#'
#' Published CCVAs usually report a sorted list of point estimates, which
#' implies a precision the data rarely support. Two complementary sources of
#' uncertainty are available:
#'
#' * `type = "bootstrap"` -- **sampling uncertainty**: units are resampled
#'   with replacement, the *whole* pipeline (component scoring, constraints,
#'   combination model) is refitted on each resample, and the original units
#'   are scored with the refitted model. For the latent model each unit's
#'   score is a prediction of an unobserved quantity, so (with
#'   `score_error = TRUE`) its factor-score standard error is added to every
#'   replicate; without it, intervals reflect parameter uncertainty only and
#'   are far too narrow.
#' * `type = "weights"` -- **weighting uncertainty** for the weighted index
#'   (a partial-identification view): component weights are drawn uniformly
#'   from the simplex, so the intervals show how much a ranking depends on
#'   the arbitrary choice of weights.
#'
#' Units are often not independent: related species share traits, and
#' species living in the same lake or river basin share their range and
#' hence their exposure. Resampling such units one at a time treats them as
#' independent replicates and makes the intervals too narrow. With
#' `cluster`, the bootstrap resamples whole groups (a cluster bootstrap):
#' every unit of a drawn group enters the replicate, so dependence within
#' groups is preserved.
#'
#' For each unit the result gives an interval for its score and for its
#' rank; for each pair of units, the probability that one is more vulnerable
#' than the other. A ranking of two units is called *robust* when that
#' probability is above `1 - robust` or below `robust`.
#'
#' @param fit A `ccva_fit`.
#' @param type `"bootstrap"` or `"weights"`.
#' @param R Number of replicates.
#' @param level Interval coverage.
#' @param robust Threshold for calling a pairwise ranking robust.
#' @param score_error For the latent model, add each unit's factor-score
#'   prediction error to the bootstrap draws (recommended).
#' @param seed Optional random seed.
#' @param progress Show a progress bar.
#' @param cluster For `type = "bootstrap"`: optional grouping of the units
#'   (for example family, genus or river basin), given as the name of a
#'   column of the data or as a vector with one value per unit. Groups are
#'   resampled with replacement instead of units. Units with a missing group
#'   form groups of their own.
#' @return An object of class `ccva_uncertainty` with `summary` (one row
#'   per unit: estimate, interval, rank interval), `pairwise` (matrix of
#'   P(row unit more vulnerable than column unit)), `robust_share`,
#'   `overlap_share` and the raw `draws`.
#' @export
#' @examples
#' fit <- ccva_combine(ccva_sim_data(), "weighted")
#' u <- ccva_uncertainty(fit, "weights", R = 500, seed = 1)
#' u
#' \donttest{
#' ub <- ccva_uncertainty(fit, "bootstrap", R = 50, seed = 1)
#' head(ub$summary)
#'
#' # cluster bootstrap: resample groups of related units
#' grp <- rep(1:60, each = 5)
#' uc <- ccva_uncertainty(fit, "bootstrap", R = 50, seed = 1, cluster = grp)
#' uc
#' }
ccva_uncertainty <- function(fit, type = c("bootstrap", "weights"),
                             R = 200, level = 0.95, robust = 0.05,
                             score_error = TRUE, seed = NULL,
                             progress = interactive(), cluster = NULL) {
  type <- match.arg(type)
  if (!inherits(fit, "ccva_fit")) stop("`fit` must be a ccva_fit.",
                                       call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  x <- fit$data
  n <- length(x$id)
  failures <- 0L
  groups <- .cluster_groups(x, cluster)
  if (!is.null(groups) && type == "weights")
    stop("`cluster` applies to type = \"bootstrap\" only.", call. = FALSE)

  if (type == "weights") {
    if (fit$method != "weighted")
      stop("type = \"weights\" applies to the weighted index only.",
           call. = FALSE)
    C <- as.matrix(fit$component_scores)
    s <- fit$details$signs[colnames(C)]
    k <- ncol(C)
    W <- matrix(stats::rexp(R * k), R, k)
    W <- W / rowSums(W)
    draws <- t(apply(W, 1, function(w)
      (C %*% (w * s) + sum(w[s < 0])) / sum(w)))
  } else {
    if (progress) pb <- utils::txtProgressBar(0, R, style = 3)
    op <- options(ccva.light_refit = TRUE)
    on.exit(options(op), add = TRUE)
    draws <- matrix(NA_real_, R, n)
    for (b in seq_len(R)) {
      rows <- if (is.null(groups)) sample.int(n, n, replace = TRUE) else
        .cluster_rows(groups)
      xb <- .subset_ccva(x, rows)
      p <- tryCatch(suppressMessages(suppressWarnings({
        fb <- do.call(ccva_combine, c(list(x = xb), fit$args))
        raw <- fb$predict_raw(x$data)
        se <- fb$details$score_se
        if (score_error && !is.null(se)) raw <- raw + stats::rnorm(n, 0, se)
        .rescale01(raw, fb$raw_range[1], fb$raw_range[2])
      })), error = function(e) NULL)
      if (is.null(p) || any(!is.finite(p))) {
        failures <- failures + 1L
      } else {
        draws[b, ] <- p
      }
      if (progress) utils::setTxtProgressBar(pb, b)
    }
    if (progress) close(pb)
    draws <- draws[stats::complete.cases(draws), , drop = FALSE]
    if (nrow(draws) < 10)
      stop("Too few successful bootstrap replicates (", nrow(draws), ").",
           call. = FALSE)
    if (failures)
      message(failures, " of ", R, " bootstrap replicates failed and were ",
              "dropped.")
  }
  colnames(draws) <- x$id

  a <- (1 - level) / 2
  ranks <- t(apply(-draws, 1, rank, ties.method = "min"))
  summ <- data.frame(
    unit = x$id, estimate = as.numeric(fit$score),
    lower = apply(draws, 2, stats::quantile, a, names = FALSE),
    upper = apply(draws, 2, stats::quantile, 1 - a, names = FALSE),
    rank = rank(-fit$score, ties.method = "min"),
    rank_low = apply(ranks, 2, stats::quantile, a, names = FALSE),
    rank_high = apply(ranks, 2, stats::quantile, 1 - a, names = FALSE),
    row.names = NULL)
  summ$width <- summ$upper - summ$lower

  # pairwise probabilities P(i more vulnerable than j)
  P <- matrix(0, n, n, dimnames = list(x$id, x$id))
  for (r in seq_len(nrow(draws))) P <- P + outer(draws[r, ], draws[r, ], ">")
  P <- P / nrow(draws)
  ut <- upper.tri(P)
  robust_share <- mean(P[ut] > 1 - robust | P[ut] < robust)
  ov <- outer(summ$lower, summ$upper, "<") & outer(summ$upper, summ$lower, ">")
  overlap_share <- mean(ov[ut])

  structure(list(summary = summ[order(summ$rank), ], pairwise = P,
                 robust_share = robust_share, overlap_share = overlap_share,
                 draws = draws, type = type, level = level, robust = robust,
                 R = nrow(draws), failures = failures, method = fit$method,
                 cluster = if (is.null(groups)) NULL else
                   list(n_groups = length(unique(groups)),
                        name = if (is.character(cluster) && length(cluster) == 1)
                          cluster else "user-supplied")),
            class = "ccva_uncertainty")
}

#' @export
print.ccva_uncertainty <- function(x, ...) {
  lab <- if (x$type == "bootstrap") {
    if (is.null(x$cluster)) "bootstrap (sampling)" else
      sprintf("cluster bootstrap over %d groups (%s)", x$cluster$n_groups,
              x$cluster$name)
  } else "random weights (weighting)"
  cat(sprintf("<ccva_uncertainty> %s, %d replicates, %d%% intervals\n",
              lab, x$R, round(100 * x$level)))
  s <- x$summary
  cat(sprintf("  mean interval width: %.3f (score range %.3f)\n",
              mean(s$width), diff(range(s$estimate))))
  cat(sprintf("  pairs with overlapping intervals: %.1f%%\n",
              100 * x$overlap_share))
  cat(sprintf("  robust pairwise rankings (P > %.2f): %.1f%%\n",
              1 - x$robust, 100 * x$robust_share))
  top <- utils::head(s, 5)
  cat("  most vulnerable:\n")
  for (i in seq_len(nrow(top)))
    cat(sprintf("    %-16s %.3f [%.3f, %.3f]  rank %d (%g-%g)\n",
                top$unit[i], top$estimate[i], top$lower[i], top$upper[i],
                top$rank[i], round(top$rank_low[i]), round(top$rank_high[i])))
  invisible(x)
}

# Cluster membership of every unit (NULL for the ordinary bootstrap)
.cluster_groups <- function(x, cluster) {
  if (is.null(cluster)) return(NULL)
  n <- length(x$id)
  g <- if (is.character(cluster) && length(cluster) == 1 && n != 1) {
    if (!cluster %in% names(x$data))
      stop("Cluster column '", cluster, "' not found in the data.",
           call. = FALSE)
    x$data[[cluster]]
  } else cluster
  if (length(g) != n)
    stop("`cluster` must name a column or give one value per unit.",
         call. = FALSE)
  g <- as.character(g)
  miss <- is.na(g)
  g[miss] <- paste0(".unit_", which(miss))
  if (length(unique(g)) < 5)
    stop("A cluster bootstrap needs at least 5 groups.", call. = FALSE)
  g
}

# Rows of one cluster-bootstrap replicate: draw as many groups as there are,
# with replacement, and take all units of each drawn group
.cluster_rows <- function(groups) {
  idx <- split(seq_along(groups), groups)
  drawn <- sample.int(length(idx), length(idx), replace = TRUE)
  unlist(idx[drawn], use.names = FALSE)
}
