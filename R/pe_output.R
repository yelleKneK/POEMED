# Assemble the tidy result table shared by all three pe_mediation_*() workers,
# so the continuous, binary, and count outcomes return one consistent schema.
# The table reports, side by side, the benchmark Wald test on the total
# indirect effect (labeled "hdmm", after Guo et al.) and the power-enhanced
# test ("pe"), as the article's data-analysis tables
# do, plus the estimated total indirect effect. No confidence interval is
# reported: the article reports none, and the Wald interval built on the
# penalized estimate under-covers in finite samples (measured 2026-09-28).
# The selected mediator set (the candidates the screen
# selected), the per-mediator table, and any per-method
# selections are variable length and not numeric, so they travel as
# attributes (keeping the value column numeric) and are surfaced through
# the print footer and the accessors pe_mediators() / pe_selection(). Not
# exported.
#
# Arguments are the already-computed scalars from a worker:
#   stat_hdmm, pval_hdmm   the Wald statistic S_n and its chi-square_q p-value
#   stat_pe, j_pe, pval_pe the PE statistic M_PE = S_n + J_m, J_m, and p-value
#   beta_hat               the estimated total indirect effect (length q)
#   selected_mediators     indices (into 1..p) of the mediators the screen
#                          selected (its estimate of the active set)
#   n, p, q                sample size, number of candidate mediators, exposures
#   outcome_label          human-readable outcome model name for the footer
#   method, error_level    the multiplicity method and its target FWER/FDR level
#   mediator_table         per-mediator data.frame (attribute), or NULL
#   selection_by_method    named list of selected sets per method (attribute), or NULL
#   pe_by_method           data.frame of the PE statistic/p-value per method
#                          (attribute), or NULL; columns method, stat_pe, j_pe,
#                          pval_pe, n_selected
#' @keywords internal
#' @noRd
.pe_assemble_output <- function(stat_hdmm, pval_hdmm, stat_pe, j_pe, pval_pe,
                                beta_hat, selected_mediators, n, p, q,
                                outcome_label, method, error_level,
                                mediator_table = NULL,
                                selection_by_method = NULL,
                                pe_by_method = NULL, tuning = NULL,
                                empty_fit = FALSE, exposure_names = NULL) {
  beta_hat <- as.numeric(beta_hat)
  # One total-indirect-effect row when there is a single exposure (the common
  # case), otherwise one row per exposure so the value column stays scalar,
  # suffixed with the exposure's column name when `X` has usable names and
  # with its position otherwise.
  if (length(beta_hat) == 1L) {
    beta_terms <- "total_indirect_effect"
  } else {
    sfx <- if (!is.null(exposure_names) &&
               length(exposure_names) == length(beta_hat))
      exposure_names else seq_along(beta_hat)
    beta_terms  <- paste0("total_indirect_effect_", sfx)
  }

  terms  <- c("stat_hdmm", "pval_hdmm", "stat_pe", "j_pe", "pval_pe",
              beta_terms,
              "n_selected_mediators", "df", "n_candidate_mediators",
              "n_observations")
  values <- c(stat_hdmm, pval_hdmm, stat_pe, j_pe, pval_pe, beta_hat,
              length(selected_mediators), q, p, n)

  out <- data.frame(term = terms, value = values, stringsAsFactors = FALSE)

  # p-value rows print to fixed decimals via the poemed_tbl format method.
  out <- .as_poemed_tbl(out, p_terms = c("pval_hdmm", "pval_pe"))
  attr(out, "outcome")             <- outcome_label
  attr(out, "selected_mediators")  <- selected_mediators
  attr(out, "method")              <- method
  attr(out, "error_level")         <- error_level
  attr(out, "mediator_table")      <- mediator_table
  attr(out, "selection_by_method") <- selection_by_method
  attr(out, "pe_by_method")        <- pe_by_method
  # The tuning record: the grid(s) searched and the HBIC-selected value(s),
  # so a reader can see which lambda produced the fit and whether it sat at
  # either end of the grid.
  attr(out, "tuning")              <- tuning
  attr(out, "empty_fit")           <- isTRUE(empty_fit)
  attr(out, "exposure_names")      <- exposure_names
  out
}

# The tuning record attached to every fit. `grid` and `selected` are the
# candidate grid and the HBIC-selected value for the full model;
# `selected_reduced` is the choice for the reduced model of the benchmark
# Wald test, searched on the same grid, when the family fits one (continuous
# outcome). `at_lower_end` and `at_upper_end` flag a full-model selection at
# the smallest or largest grid value, where the grid's end rather than an
# interior HBIC minimum may have determined the fit. Not exported.
#' @keywords internal
#' @noRd
.pe_tuning_record <- function(grid, selected, selected_reduced = NULL) {
  rec <- list(lambda_grid = grid, lambda_selected = selected,
              at_lower_end = isTRUE(selected <= min(grid)),
              at_upper_end = isTRUE(selected >= max(grid)))
  if (!is.null(selected_reduced)) rec$lambda_selected_reduced <- selected_reduced
  rec
}

# The warning raised when the penalized fit that HBIC selects from the grid
# contains no mediator. The condition tested is only that the selected fit
# is empty: smaller grid values may well have kept mediators and lost to the
# empty fit on HBIC, so the message names the selected lambda and does not
# claim that every grid value was empty. Under a complete null an empty fit
# is the expected result. `selected` is the HBIC-selected lambda (NULL
# leaves it out of the message). The warning carries the class
# "poemed_empty_fit" so the Monte Carlo functions can count these fits
# without printing hundreds of warnings. Not exported.
#' @keywords internal
#' @noRd
.pe_empty_fit_warning <- function(grid, selected = NULL) {
  at <- if (is.null(selected)) "" else sprintf("lambda = %.3g, ", selected)
  msg <- sprintf(paste0(
    "The penalized fit that HBIC selected (%sfrom %d values of `lambda_grid`, ",
    "%.3g to %.3g) contains no mediator, so the benchmark and power-enhanced ",
    "tests were not computed; both are reported as 0 with p-values of 1, ",
    "POEMED's reporting convention for an empty selection. An empty fit is ",
    "expected when no mediator is active, but a p-value of 1 here is not ",
    "evidence for the null; see the `lambda_grid` argument of ?pe_mediation."),
    at, length(grid), min(grid), max(grid))
  warning(structure(class = c("poemed_empty_fit", "warning", "condition"),
                    list(message = msg, call = NULL)))
}

# The warning raised when the HBIC-selected lambda is the smallest or largest
# value of the grid searched (for a continuous outcome the reduced model's
# choice is checked too). HBIC is minimized over the grid alone, so a minimum at an
# end may lie beyond it, or between the end and its neighbor when the grid is
# coarse; the remedy is to extend the grid past that end and to space its
# values more finely. Called only for a non-empty fit; the empty-fit warning
# covers the other case. The warning carries the class "poemed_grid_boundary"
# so the Monte Carlo functions and WHO_mediation_analysis(), which search a
# fixed grid by design, can muffle it. Not exported.
#' @keywords internal
#' @noRd
.pe_grid_boundary_warning <- function(tuning) {
  g <- tuning$lambda_grid
  end_of <- function(selected) {
    if (isTRUE(selected <= min(g))) "smallest"
    else if (isTRUE(selected >= max(g))) "largest" else NULL
  }
  full <- end_of(tuning$lambda_selected)
  has_reduced <- !is.null(tuning$lambda_selected_reduced)
  red <- if (has_reduced) end_of(tuning$lambda_selected_reduced) else NULL
  if (is.null(full) && is.null(red)) return(invisible(FALSE))
  phrase <- function(selected, end, model = NULL)
    paste0(sprintf("lambda = %.3g, the %s value of `lambda_grid`", selected, end),
           if (!is.null(model)) paste0(", for ", model))
  hits <- if (!has_reduced) phrase(tuning$lambda_selected, full)
  else if (!is.null(full) && identical(full, red) &&
           isTRUE(all.equal(tuning$lambda_selected, tuning$lambda_selected_reduced)))
    phrase(tuning$lambda_selected, full,
           "the full model and the reduced model of the benchmark test")
  else c(if (!is.null(full)) phrase(tuning$lambda_selected, full, "the full model"),
         if (!is.null(red)) phrase(tuning$lambda_selected_reduced, red,
                                   "the reduced model of the benchmark test"))
  msg <- paste0(
    "HBIC selected ", paste(hits, collapse = " and "), ". The criterion is ",
    "minimized over the grid alone, so its minimum may lie beyond that end or ",
    "between the end and its neighbor. Consider extending the grid past it and ",
    "using a finer partition (more values), and compare the selected mediators ",
    "and p-values across grids.")
  warning(structure(class = c("poemed_grid_boundary", "warning", "condition"),
                    list(message = msg, call = NULL)))
  invisible(TRUE)
}

# Post-inference finalization shared by the three workers. Given the benchmark
# Wald pieces (statistic, p-value, beta, and its asymptotic covariance) and the
# selected mediators, it (1) runs the power-enhancement screen for each
# requested multiplicity method (the first is primary; the rest, if any, are
# recorded for comparison), and (2) assembles the tidy output. Centralizing
# this keeps the workers identical past the family-specific fit. Not exported.
#
# The covariance var_beta is the asymptotic covariance of sqrt(n)*beta_hat (the
# matrix the Wald statistic S_n = n beta' var_beta^{-1} beta inverts), so the
# variance of beta_hat itself is var_beta / n. It is passed through but not
# reported: the Wald interval it would give was found to under-cover in finite
# samples for binary and count outcomes (the penalized estimate is shrunk), so
# the package reports no interval for the total indirect effect.
#' @keywords internal
#' @noRd
.pe_finalize <- function(X, Y, M_A, S, A, stat_hdmm, pval_hdmm, beta_hat,
                         var_beta, p, n, q, y_family, methods, error_level,
                         outcome_label, tuning = NULL,
                         exposure_names = NULL) {
  # Screen with each requested method off the SAME penalized fit (the penalized
  # estimation, the expensive part, is already done; the screen only refits the
  # cheap path regressions). The first method drives the primary output.
  pe_list <- lapply(methods, function(m)
    .pe_component(X, Y, M_A, S, A, Sn = stat_hdmm, p_total = p, n = n, q = q,
                  y_family = y_family, method = m, error_level = error_level))
  names(pe_list) <- methods
  primary <- pe_list[[1L]]
  selection_by_method <- if (length(methods) > 1L)
    lapply(pe_list, function(pc) pc$selected_mediators) else NULL

  # The PE statistic, J_m, and p-value differ by multiplicity method (each
  # method screens a different selected set into J_m), so record all of them.
  # This is what pe_selection() reports side by side and the article's
  # "full" data-analysis tables show as the PE / PE_BH / PE_BY columns.
  pe_by_method <- data.frame(
    method   = methods,
    stat_pe  = vapply(pe_list, function(pc) pc$stat_pe, numeric(1)),
    j_pe     = vapply(pe_list, function(pc) pc$j_pe, numeric(1)),
    pval_pe  = vapply(pe_list, function(pc) pc$pval_pe, numeric(1)),
    n_selected = vapply(pe_list, function(pc) length(pc$selected_mediators),
                        integer(1)),
    row.names = NULL, stringsAsFactors = FALSE)

  .pe_assemble_output(stat_hdmm, pval_hdmm, primary$stat_pe, primary$j_pe,
                      primary$pval_pe, beta_hat,
                      primary$selected_mediators, n, p, q, outcome_label,
                      methods[1L], error_level,
                      mediator_table = primary$mediator_table,
                      selection_by_method = selection_by_method,
                      pe_by_method = pe_by_method, tuning = tuning,
                      exposure_names = exposure_names)
}

# Map mediator indices in a finished result from positions in the (possibly
# reduced) mediator matrix back to positions in the user's original `M`, and
# attach the mediator names when `M` had them. `kept` is the vector of
# surviving original-column indices when drop_constant dropped zero-variance
# columns (a reduced index j becomes kept[j]); NULL when nothing was dropped.
# `mediator_names` are the column names of the user's original `M` (NULL when
# unnamed); with names, the mediator table gains a `name` column and the
# names travel in the "mediator_names" attribute for the print footer,
# pe_selection(), and plot(). Not exported.
#' @keywords internal
#' @noRd
.pe_relabel_mediators <- function(out, kept, mediator_names = NULL) {
  if (!is.null(kept)) {
    am <- attr(out, "selected_mediators")
    if (length(am)) attr(out, "selected_mediators") <- kept[am]
    sbm <- attr(out, "selection_by_method")
    if (!is.null(sbm))
      attr(out, "selection_by_method") <-
        lapply(sbm, function(idx) if (length(idx)) kept[idx] else idx)
    mt <- attr(out, "mediator_table")
    if (!is.null(mt) && nrow(mt)) {
      mt$mediator <- kept[mt$mediator]
      attr(out, "mediator_table") <- mt
    }
  }
  attr(out, "mediator_names") <- mediator_names
  mt <- attr(out, "mediator_table")
  if (!is.null(mediator_names) && !is.null(mt)) {
    rest <- setdiff(names(mt), "mediator")
    mt <- cbind(mt[, "mediator", drop = FALSE],
                data.frame(name = mediator_names[mt$mediator],
                           stringsAsFactors = FALSE),
                mt[, rest, drop = FALSE])
    rownames(mt) <- NULL
    attr(out, "mediator_table") <- mt
  }
  out
}

# The result returned when the HBIC-selected penalized fit contains no
# candidate-active mediator (the selected set S-hat is empty). With no
# selected mediator there is no indirect effect to test: J_m = 0 (the article
# states only that J_m = 0 when the selection is empty), and the benchmark
# statistic is reported as 0 with a p-value of 1, POEMED's
# reporting convention (an empty selection is a non-rejection). The result
# keeps the full row schema, carries `empty_fit = TRUE`, and is announced by
# a classed warning so it cannot be mistaken for evidence for the null.
#
# `methods` is the vector of multiplicity methods requested (length > 1 when
# report_all_methods = TRUE); with no selected mediators every method gives the
# same all-zero result, so they share J_m = 0, stat_pe = 0, pval_pe = 1.
#' @keywords internal
#' @noRd
.pe_empty_output <- function(n, p, q, outcome_label, method, error_level,
                             methods = method,
                             tuning = NULL, exposure_names = NULL) {
  pe_by_method <- data.frame(
    method = methods, stat_pe = 0, j_pe = 0, pval_pe = 1, n_selected = 0L,
    row.names = NULL, stringsAsFactors = FALSE)
  selection_by_method <- if (length(methods) > 1L)
    stats::setNames(rep(list(integer(0)), length(methods)), methods) else NULL
  if (!is.null(tuning))
    .pe_empty_fit_warning(tuning$lambda_grid, tuning$lambda_selected)
  .pe_assemble_output(
    stat_hdmm = 0, pval_hdmm = 1, stat_pe = 0, j_pe = 0, pval_pe = 1,
    beta_hat = rep(0, q),
    selected_mediators = integer(0), n = n, p = p, q = q,
    outcome_label = outcome_label, method = methods[1L],
    error_level = error_level,
    mediator_table = data.frame(
      mediator = integer(0), t_outcome = numeric(0), t_exposure = numeric(0),
      screen_p = numeric(0), selected = logical(0)),
    selection_by_method = selection_by_method,
    pe_by_method = pe_by_method, tuning = tuning, empty_fit = TRUE,
    exposure_names = exposure_names)
}
