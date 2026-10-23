## Rule-based combination (Foden et al. 2013) --------------------------------
##
## Each indicator flags a unit as "high" (exposure, sensitivity) or "low"
## (adaptive capacity). Binary or three-level indicators coded 0 / (0.5) / 1
## (reversed when their direction is -1) flag at 1 (exposure, sensitivity) or 0 (adaptive
## capacity); continuous indicators flag in the top `threshold` share
## (exposure, sensitivity) or bottom `threshold` share (adaptive capacity).
## A component is flagged when any of its indicators is. A unit is highly
## vulnerable when all components are flagged.

.rule_levels <- c("highly vulnerable", "high latent risk",
                  "potential adapter", "potential persister", "lower risk")

.rule_flag_label <- function(comp) {
  lab <- c(exposure = "high exposure", sensitivity = "high sensitivity",
           adaptive = "low adaptive capacity")
  unname(lab[comp])
}

.is_coded <- function(v) {
  u <- unique(v[!is.na(v)])
  length(u) > 0 && all(u %in% c(0, 0.5, 1))
}

.comb_rule <- function(x, comps, cons, threshold, unknown) {
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0 ||
      threshold >= 1)
    stop("`threshold` must be a single number between 0 and 1.",
         call. = FALSE)
  notes <- character(0)
  if (length(cons$comp_bounds) || length(cons$equal) ||
      length(cons$ind_bounds) || !is.null(cons$moderator))
    notes <- "bounds, equality and moderator rules are not expressible in the rule-based method and were not applied"

  used <- comps
  for (cmp in names(cons$comp_fixed))
    if (cons$comp_fixed[[cmp]] == 0) used <- setdiff(used, cmp)
  dropped <- names(Filter(function(v) v == 0, cons$ind_fixed))

  # per-indicator rules, learned from the data
  rules <- list()
  for (cmp in used) {
    vars <- setdiff(x$indicators[[cmp]], dropped)
    X <- .oriented(x, cmp)[, vars, drop = FALSE]
    low_is_risk <- cmp == "adaptive"
    for (v in vars) {
      xv <- X[, v]
      if (.is_coded(x$data[[v]])) {
        rules[[v]] <- list(component = cmp, type = "coded",
                           cut = if (low_is_risk) 0 else 1,
                           low = low_is_risk)
      } else {
        q <- stats::quantile(xv, if (low_is_risk) threshold else 1 - threshold,
                             names = FALSE, na.rm = TRUE)
        rules[[v]] <- list(component = cmp, type = "quantile", cut = q,
                           low = low_is_risk)
      }
    }
  }

  flag_matrix <- function(data) {
    sapply(names(rules), function(v) {
      r <- rules[[v]]
      xv <- if (r$type == "coded") {
        # coded 0 / 0.5 / 1 on the raw scale; reverse when direction is -1
        if (x$direction[[v]] < 0) 1 - data[[v]] else data[[v]]
      } else data[[v]] * x$direction[[v]]
      f <- if (r$low) xv <= r$cut else xv >= r$cut
      f[is.na(f)] <- unknown == "pessimistic"
      f
    })
  }
  classify <- function(data) {
    Fi <- flag_matrix(data)
    if (is.null(dim(Fi))) Fi <- matrix(Fi, nrow = 1,
                                       dimnames = list(NULL, names(rules)))
    comp_of <- vapply(rules, `[[`, character(1), "component")
    Fc <- sapply(used, function(cmp)
      apply(Fi[, comp_of == cmp, drop = FALSE], 1, any))
    if (is.null(dim(Fc))) Fc <- matrix(Fc, nrow = 1,
                                       dimnames = list(NULL, used))
    list(indicator = Fi, component = Fc)
  }

  cl <- classify(x$data)
  Fc <- cl$component; Fi <- cl$indicator
  rownames(Fc) <- x$id
  score <- rowMeans(Fc)
  has <- function(cmp) if (cmp %in% used) Fc[, cmp] else rep(FALSE, nrow(Fc))
  E <- has("exposure"); S <- has("sensitivity"); A <- has("adaptive")
  if (!"adaptive" %in% used) A <- rep(TRUE, nrow(Fc))
  catg <- ifelse(E & S & A, .rule_levels[1],
         ifelse(S & A & !E, .rule_levels[2],
         ifelse(S & E & !A, .rule_levels[3],
         ifelse(E & A & !S, .rule_levels[4], .rule_levels[5]))))
  category <- factor(catg, levels = .rule_levels)

  # drivers: how often each indicator flags, and how often it is the only
  # indicator flagging its component (the sensitivity of the result to it)
  comp_of <- vapply(rules, `[[`, character(1), "component")
  n_in_comp <- sapply(used, function(cmp)
    rowSums(Fi[, comp_of == cmp, drop = FALSE]))
  if (is.null(dim(n_in_comp))) n_in_comp <- matrix(n_in_comp, nrow = 1)
  colnames(n_in_comp) <- used
  drivers <- data.frame(
    component = comp_of, indicator = names(rules),
    type = vapply(rules, `[[`, character(1), "type"),
    cut = vapply(rules, function(r) r$cut, numeric(1)),
    n_flagged = colSums(Fi),
    n_only = vapply(names(rules), function(v)
      sum(Fi[, v] & n_in_comp[, comp_of[[v]]] == 1), numeric(1)),
    n_unknown = vapply(names(rules), function(v) sum(is.na(x$data[[v]])),
                       numeric(1)),
    row.names = NULL)

  list(score = as.numeric(score), model = list(rules = rules),
       details = list(flags = Fc, category = stats::setNames(category, x$id),
                      drivers = drivers, threshold = threshold,
                      unknown = unknown, rules = rules),
       notes = notes,
       predict = function(newdata) rowMeans(classify(newdata)$component))
}
