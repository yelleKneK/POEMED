# Fit one subgroup with the preprocessing of the analysis behind the
# article's tables: standardize the exposure and mediators, center the
# outcome, and leave the confounders (region/income indicators and year)
# UNSCALED, as that analysis does (the article's text says only that
# covariates are standardized). Returns one summary row.
# An empty penalized fit (the fit HBIC selects contains no mediator) keeps its
# convention values (p-values of 1, no selected mediator) but its classed
# warning is muffled here and the row carries a "notes" attribute instead,
# so the caller can name the group in one warning; the attribute travels
# with the row out of a forked child when cores > 1, where a warning would
# be lost. Not exported.
#' @keywords internal
#' @noRd
.WHO_fit_group <- function(grouping, group, des, lambda_grid, error_level,
                           full_table = FALSE) {
  # Some indicators are constant within a small region/income cell; drop them
  # (they carry no signal and break standardization) and keep the surviving
  # indicator codes so the reported selected mediators map back correctly. A
  # column with a missing or non-finite value has sd NA and is dropped here
  # too; WHO_mediation_analysis() names those (see .WHO_nonfinite()).
  keep <- which(apply(des$M, 2L, stats::sd) > 0)
  M <- des$M[, keep, drop = FALSE]
  ind <- des$indicators[keep]
  fit <- withCallingHandlers(
    pe_mediation_linear(
      scale(des$X), des$Y - mean(des$Y), scale(M), Z = des$Z,
      scale = FALSE, error_level = error_level,
      report_all_methods = full_table,
      lambda_grid = lambda_grid),
    poemed_empty_fit = function(w) invokeRestart("muffleWarning"),
    # The grid is the article's fixed grid: in 67 of the 128 cells the five
    # outcomes fit, the HBIC choice for the full or the reduced model (54
    # cells each) is its smallest value, a property of that analysis, not
    # something to flag cell by cell.
    poemed_grid_boundary = function(w) invokeRestart("muffleWarning"))
  pval_hdmm <- fit$value[fit$term == "pval_hdmm"]
  codes <- function(idx) if (length(idx))
    paste(ind[idx], collapse = ", ") else "none"

  if (!full_table) {
    sel_idx <- attr(fit, "selected_mediators")
    row <- data.frame(
      grouping = grouping, group = group, n_countries = des$n_countries,
      pval_hdmm = pval_hdmm,
      pval_pe   = fit$value[fit$term == "pval_pe"],
      n_selected  = length(sel_idx),
      selected_mediators = codes(sel_idx),
      stringsAsFactors = FALSE)
  } else {
    # full_table: report the PE p-value and the selected mediator codes under
    # each multiplicity method (Bonferroni / BH / BY), as the article's
    # extended data-analysis tables do.
    pbm <- attr(fit, "pe_by_method")
    sbm <- attr(fit, "selection_by_method")
    meth_key <- c(Bonferroni = "bonferroni", BH = "bh", BY = "by")
    row <- data.frame(grouping = grouping, group = group,
                      n_countries = des$n_countries,
                      pval_hdmm = pval_hdmm, stringsAsFactors = FALSE)
    for (m in names(meth_key)) {
      k <- meth_key[[m]]
      row[[paste0("pval_pe_", k)]] <-
        if (!is.null(pbm)) pbm$pval_pe[match(m, pbm$method)] else NA_real_
      idx <- if (!is.null(sbm)) sbm[[m]] else integer(0)
      row[[paste0("selected_", k)]] <- codes(idx)
    }
  }
  if (isTRUE(attr(fit, "empty_fit"))) {
    lambda <- sprintf("lambda = %.3g", attr(fit, "tuning")$lambda_selected)
    attr(row, "notes") <- .WHO_note(
      "empty",
      paste0("The penalized fit that HBIC selected (", lambda, ") contains ",
             "no mediator, so the benchmark and power-enhanced tests were ",
             "not computed; the p-values of 1 are POEMED's reporting ",
             "convention for an empty selection, not evidence for the null."),
      detail = lambda)
  }
  row
}

# The indicators with a missing or non-finite value among a design's rows.
# These are exactly the columns whose standard deviation is NA, which
# .WHO_fit_group() drops along with the constant ones; this only names them.
# Not exported.
#' @keywords internal
#' @noRd
.WHO_nonfinite <- function(des) {
  des$indicators[is.na(apply(des$M, 2L, stats::sd))]
}

# One per-group note: the stage it arose at, a self-contained message for
# the "notes" attribute, and the short detail the aggregate report lists
# after the group's name. Not exported.
#' @keywords internal
#' @noRd
.WHO_note <- function(stage, message, detail = sub("[.]$", "", message)) {
  data.frame(stage = stage, message = message, detail = detail,
             stringsAsFactors = FALSE)
}

#' Run the Article's Empirical Mediation Analysis
#'
#' Runs the article's real-data analysis for one health outcome across the
#' groupings it reports: a global model using all WHO members, then
#' separate models within each WHO region, within each World Bank income
#' group, and (optionally) within each region-by-income cell. For every
#' grouping it reports the benchmark Wald test and the power-enhanced test
#' on the total indirect effect of health expenditure, together with the
#' individual mediators the PE screen selects. The result has
#' the layout of the article's data-analysis tables (Table 1 for infant
#' mortality and Table 2 for undernourishment; Tables S.8 to S.12 of its
#' supplement give the extended layout and the other three outcomes),
#' computed from the shipped [WHO_health_mediation] data.
#'
#' The output is the package's \strong{benchmark fit}: the test suite pins
#' the selected sets and the benchmark p-values of the article's tables
#' under the function's default grid, so a change to the estimation code
#' that moved them would be caught, and the fit serves as a fixed reference
#' point for the method.
#'
#' @details
#' Each subgroup is fit with the preprocessing of the analysis behind the
#' article's tables: the exposure and the 57 mediators are standardized,
#' the outcome is centered, and the confounders (the region and income
#' indicators and the year) are left unscaled. The article's text says
#' that covariates are standardized; the unscaled confounders are a choice
#' of that analysis.
#' Within a group, an indicator that is constant there (common in the
#' small region-by-income cells) is dropped before fitting, since it
#' carries no information. An indicator with a missing or non-finite value
#' among the group's rows is dropped as well, because the method has no
#' missing-data mechanism; those are named in a message, since they change
#' the set of candidate mediators.
#'
#' Four kinds of event are reported once, after all groups are fit, and
#' recorded in the `"notes"` attribute (a `data.frame` with `group`,
#' `stage`, and `message`, or `NULL` when there are none). The `stage`
#' column says which:
#'
#' * `"design"`: a group with too few complete observations (some
#'   region-by-income cells are empty or nearly so), or whose rows lack a
#'   country code, a year, or a region or income group that varies within
#'   the group, is skipped and reported as an `NA` row, named in a message.
#' * `"missing"`: indicators with missing or non-finite values were dropped
#'   in the group before fitting, named in a message.
#' * `"empty"`: the penalized fit that HBIC selected contains no mediator.
#'   The row keeps the values of POEMED's reporting convention for an empty
#'   selection (p-values of 1, no selected mediator), which are not evidence
#'   for the null, and the group is named, with the lambda HBIC selected,
#'   in a warning of class `poemed_empty_fit` at any value of `cores`.
#' * `"fit"`: the fit failed; the group is reported as an `NA` row and named
#'   in a warning that quotes the error.
#'
#' Because every group fits a penalized model over a grid of tuning
#' parameters, a full run with all groupings fits dozens of models and
#' takes about ten seconds on a current laptop; start with the default
#' groupings, or a single outcome, before scaling up. In 67 of the 128 cells
#' that the five outcomes fit, the HBIC choice for the full penalized model
#' or for the reduced model of the benchmark Wald test (54 cells each) is
#' 0.1, the grid's smallest value; that is a property of the article's
#' analysis on its fixed grid, so the function does not raise the
#' single-fit `poemed_grid_boundary` warning for it. To see a cell's tuning
#' record, fit the cell by hand with [pe_mediation_linear()], as the
#' reproduction vignette shows.
#'
#' The tuning-parameter grid is `seq(0.1, 2, length.out = 100)` for the
#' global model and for every subgroup model, the grid behind the article's
#' Tables 1, 2, and S.8 to S.12 (Section 5 of Yu and Kelley, in press);
#' the article's text itself states only that the tuning parameter is
#' chosen by the HBIC criterion. The article's real-data analysis was run
#' with this package ([pe_mediation_linear()] at that grid, with the
#' preprocessing described above); `WHO_mediation_analysis()` automates
#' those calls. Under this default the function reproduces those tables. The selected set agrees with the printed one in every
#' fitted cell under Bonferroni, BY, and BH, with one exception (the BH set
#' for under-five mortality in the High income group, where the article's
#' set also holds `gge_gdp`). The benchmark p-value agrees to four decimals
#' in all but a handful of cells: infant mortality in SEAR (0.6483 here,
#' 0.6482 printed) and in AMR-High (0.0726 here, 0.0725 printed) differ by
#' rounding in the fourth decimal, and a few p-values at or below 0.0001,
#' which the article prints to one significant figure, come out somewhat
#' smaller here. Exact p-values depend on the grid and on the unscaled
#' confounders. Six of the region-by-income cells that the article's
#' tables fill are not fit at all and are reported as `NA` (the cells the
#' tables list with 0 countries are `NA` too, skipped before fitting):
#'
#' * four cells in which two indicators are exact linear partners and the
#'   fit keeps both, so a matrix it inverts is singular: infant mortality
#'   and under-five mortality in EMR-Lower-middle (`dom_che` and `ext_che`
#'   sum to 100 there), and life expectancy and under-five mortality in
#'   AFR-Upper-middle (`ext_che` with `dom_che`, and `shi_che` with
#'   `chi_che`). The article's values are reproduced by dropping one
#'   partner; see the recipe below.
#' * two cells of the undernourishment outcome, EUR-High and WPR-High, in
#'   which the outcome is 2.5 in every country-year (the lower reporting
#'   bound of the source), so there is no variation to mediate. The
#'   function reports the fit as failed with the message that `Y` is
#'   constant; the article prints a p-value of 1 there.
#'
#' Many of the 57 indicators are exact or near-exact linear combinations of
#' one another, and the penalized selection chooses among such partners.
#' As a result the benchmark p-values, and occasionally the selected
#' mediators, can change when the indicator columns of `data` are
#' reordered, even between two indicators that are exactly collinear within
#' a group. The results reported are those for the column order of the
#' shipped data.
#'
#' A group whose fit fails with a singular-matrix error usually holds two
#' indicators that are exact linear partners within that group (for
#' example, `dom_che` and `ext_che` sum to 100 in the Eastern Mediterranean
#' lower-middle-income group; `stats::cor()` of 1 or -1 among the columns
#' of the group's `M` finds such pairs). Such a group can be fit by hand
#' after dropping one member of the pair. Which member is dropped matters,
#' since the two fits select among different columns (dropping `ext_che`
#' instead in the infant-mortality cell below gives 0.0786). The article's
#' analysis removed `dom_che` and `vfa_che` in the two EMR-Lower-middle
#' cells (`cfa_che` and `vfa_che` are near-exact partners there as well),
#' giving benchmark p-values of 0.0837 for infant mortality and 0.0914 for
#' under-five mortality with no mediator selected, and `shi_che` and
#' `ext_che` in the two AFR-Upper-middle cells, giving 0.1311 for life
#' expectancy with no mediator selected and 0.0389 for under-five
#' mortality with `pvtd_ncu2021_pc` selected. With the default grid, the
#' infant-mortality cell is:
#' \preformatted{
#'   cell <- WHO_mediation_design("imr", region = "EMR",
#'                                income = "Lower-middle")
#'   keep <- apply(cell$M, 2, sd) > 0 &
#'     !(cell$indicators %in% c("dom_che", "vfa_che"))
#'   grid <- seq(0.1, 2, length.out = 100)
#'   pe_mediation(scale(cell$X), cell$Y - mean(cell$Y),
#'                scale(cell$M[, keep]), Z = cell$Z,
#'                outcome = "continuous", scale = FALSE,
#'                lambda_grid = grid)
#' }
#' which prints the benchmark p-value 0.0837 of the article's Table 1.
#'
#' @param outcome Which health outcome to analyze: `"imr"`, `"u5mr"`,
#'   `"leb"`, `"lbw"`, or `"pou"`. See [WHO_health_mediation].
#' @param groupings Character vector of groupings to run, any of
#'   `"global"`, `"region"`, `"income"`, and `"region_income"`. Default
#'   `c("global", "region", "income")` (the region-by-income cells are the
#'   slowest and are opt-in).
#' @param error_level Target familywise error rate for the
#'   mediator-screening step. Default 0.05, as in Yu and Kelley (in press).
#' @param lambda_grid Tuning-parameter grid for the subgroup penalized fits
#'   (regional, income, and region-by-income models). Default
#'   `seq(0.1, 2, length.out = 100)`, the grid behind the article's
#'   real-data tables (Tables 1, 2, and S.8 to S.12 of Yu and Kelley, in
#'   press); it differs from the default of [pe_mediation()] because the
#'   function reproduces the article's analysis as run.
#' @param lambda_grid_global Tuning-parameter grid for the global
#'   (all-country) model only. Defaults to `lambda_grid`: the article's
#'   tables use the same grid for the global model and the subgroup
#'   models.
#' @param data The source data, in the shape of [WHO_health_mediation]; see
#'   [WHO_mediation_design()] for the columns it must carry. Defaults to
#'   [WHO_health_mediation].
#' @param cores Number of CPU cores; values above 1 fit the subgroups in
#'   parallel via the base \pkg{parallel} package (Unix only). Default 1.
#' @param verbose Logical; if `TRUE`, print each group as it is fit.
#'   Default `FALSE`.
#' @param full_table Logical; if `TRUE`, return the article's
#'   \emph{extended} table layout, with the power-enhanced p-value and the
#'   selected mediator codes reported separately under each multiplicity
#'   method (Bonferroni for FWER, BH and BY for FDR) instead of a single
#'   primary method. Default `FALSE`.
#'
#' @return A `data.frame` with one row per fitted subgroup. With
#'   `full_table = FALSE` (the default) the columns are `grouping`
#'   (global / region / income / region_income), `group` (the specific
#'   group, for example `AFR` or `Low`), `n_countries` (number of countries;
#'   `NA` for a group skipped at the design stage: one with fewer than 10
#'   complete observations, which includes every group the article's tables
#'   list with 0 countries, or one whose rows lack a key the group needs),
#'   `pval_hdmm` and `pval_pe` (the benchmark and power-enhanced p-values
#'   for the total indirect effect), `n_selected` (number of mediators the
#'   screen selected), and `selected_mediators` (their codes, or `"none"`). With
#'   `full_table = TRUE` the single `pval_pe`/`n_selected`/`selected_mediators`
#'   columns are replaced by `pval_pe_bonferroni`, `pval_pe_bh`,
#'   `pval_pe_by` and the matching `selected_bonferroni`, `selected_bh`,
#'   `selected_by` code columns. The result is a `poemed_tbl` with the
#'   attributes `"outcome"`, `"full_table"`, and `"notes"` (see Details);
#'   call [summary()][summary.poemed_tbl] on it for a cross-group digest of
#'   where the power-enhanced test detects mediation the benchmark misses.
#'
#' @inheritSection POEMED-package How to Cite
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @references
#' Yu, X., & Kelley, K. (in press). Power Enhancement in
#' High-Dimensional Heterogeneous Mediation Analysis. \emph{Journal of the
#' American Statistical Association}.
#'
#' @seealso [WHO_mediation_design()] for a single design,
#'   [WHO_health_mediation] for the data, [pe_mediation()] for the test.
#'
#' @family WHO mediation data
#'
#' @examples
#' # The global infant-mortality model, the first row of the article's
#' # Table 1 (under a second): benchmark p-value 0.0644, power-enhanced
#' # p-value below 0.0001, and gge_gdp selected. Add "region" and
#' # "income" to groupings for the rest of the table.
#' WHO_mediation_analysis("imr", groupings = "global")
#'
#' @export
WHO_mediation_analysis <- function(outcome = c("imr", "u5mr", "leb",
                                               "lbw", "pou"),
                                   groupings = c("global", "region", "income"),
                                   error_level = 0.05,
                                   lambda_grid = seq(0.1, 2, length.out = 100),
                                   lambda_grid_global = lambda_grid,
                                   data = WHO_health_mediation,
                                   cores = 1L, verbose = FALSE,
                                   full_table = FALSE) {
  outcome <- .WHO_match_arg(outcome, c("imr", "u5mr", "leb", "lbw", "pou"),
                            "outcome")
  groupings <- .WHO_match_arg(groupings,
                              c("global", "region", "income", "region_income"),
                              "groupings", several.ok = TRUE)
  .validate_error_level(error_level)
  .validate_lambda_grid(lambda_grid, "lambda_grid")
  .validate_lambda_grid(lambda_grid_global, "lambda_grid_global")
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
      cores < 1 || cores != round(cores))
    stop("`cores` must be a single positive whole number.", call. = FALSE)
  for (nm in c("verbose", "full_table")) {
    v <- get(nm)
    if (!is.logical(v) || length(v) != 1L || is.na(v))
      stop(sprintf("`%s` must be TRUE or FALSE.", nm), call. = FALSE)
  }
  .WHO_validate_data(data, outcome)

  regions <- c("AFR", "AMR", "EMR", "EUR", "SEAR", "WPR")
  incomes <- c("Low", "Lower-middle", "Upper-middle", "High")

  # Build the (grouping, group, region-filter, income-filter) work list from
  # the requested groupings.
  jobs <- list()
  if ("global" %in% groupings)
    jobs <- c(jobs, list(list(g = "global", lab = "ALL", reg = NULL, inc = NULL)))
  if ("region" %in% groupings)
    jobs <- c(jobs, lapply(regions,
              function(r) list(g = "region", lab = r, reg = r, inc = NULL)))
  if ("income" %in% groupings)
    jobs <- c(jobs, lapply(incomes,
              function(i) list(g = "income", lab = i, reg = NULL, inc = i)))
  if ("region_income" %in% groupings)
    jobs <- c(jobs, unlist(lapply(regions, function(r) lapply(incomes,
              function(i) list(g = "region_income", lab = paste(r, i, sep = "-"),
                               reg = r, inc = i))), recursive = FALSE))

  # Fit one grouping cell, returning a one-row summary. Cells with too few
  # observations (WHO_mediation_design() errors) or a failed fit are recorded as
  # NA rather than aborting the whole run; indicators dropped for missing
  # values and an empty penalized fit leave the row's values as they are.
  # Each event travels as a "notes" attribute on the row (so it survives the
  # fork when cores > 1) and is reported once, after the loop.
  fit_job <- function(job) {
    if (verbose) message("Fitting ", outcome, " / ", job$lab, " ...")
    na_row <- function(n_countries) {
      if (!full_table)
        return(data.frame(
          grouping = job$g, group = job$lab, n_countries = n_countries,
          pval_hdmm = NA_real_, pval_pe = NA_real_, n_selected = NA_integer_,
          selected_mediators = NA_character_, stringsAsFactors = FALSE))
      data.frame(
        grouping = job$g, group = job$lab, n_countries = n_countries,
        pval_hdmm = NA_real_,
        pval_pe_bonferroni = NA_real_, selected_bonferroni = NA_character_,
        pval_pe_bh = NA_real_, selected_bh = NA_character_,
        pval_pe_by = NA_real_, selected_by = NA_character_,
        stringsAsFactors = FALSE)
    }
    des <- tryCatch(
      WHO_mediation_design(outcome, region = job$reg, income = job$inc,
                           data = data),
      error = function(e) e)
    if (inherits(des, "error")) {
      row <- na_row(NA_integer_)
      attr(row, "notes") <- .WHO_note("design", conditionMessage(des))
      return(row)
    }
    dropped <- .WHO_nonfinite(des)
    missing_note <- if (length(dropped))
      .WHO_note("missing",
                paste0("Dropped before fitting for missing or non-finite ",
                       "values: ", paste(dropped, collapse = ", "), "."),
                detail = paste(dropped, collapse = ", "))
    grid <- if (job$g == "global") lambda_grid_global else lambda_grid
    res <- tryCatch(
      .WHO_fit_group(job$g, job$lab, des, grid, error_level,
                     full_table = full_table),
      error = function(e) e)
    if (inherits(res, "error")) {
      row <- na_row(des$n_countries)
      attr(row, "notes") <- rbind(missing_note,
                                  .WHO_note("fit", conditionMessage(res)))
      return(row)
    }
    attr(res, "notes") <- rbind(missing_note, attr(res, "notes"))
    res
  }
  rows <- .pe_lapply(jobs, fit_job, cores = cores)

  # Collect the per-cell notes (skipped designs, dropped indicators, empty
  # fits, failed fits) and report them once: skips and drops as messages,
  # empty and failed fits as warnings.
  notes <- do.call(rbind, lapply(seq_along(rows), function(i) {
    nt <- attr(rows[[i]], "notes")
    if (is.null(nt) || !nrow(nt)) return(NULL)
    data.frame(group = jobs[[i]]$lab, nt, stringsAsFactors = FALSE)
  }))
  .WHO_report_notes(notes)
  if (!is.null(notes)) {
    notes$detail <- NULL
    rownames(notes) <- NULL
  }

  out <- do.call(rbind, lapply(rows, function(r) { attr(r, "notes") <- NULL; r }))
  rownames(out) <- NULL
  # Return the package's general poemed_tbl so summary() gives the cross-group
  # benchmark-vs-PE digest; it still behaves as an ordinary data.frame.
  out <- .as_poemed_tbl(out)
  attr(out, "outcome") <- outcome
  attr(out, "full_table") <- full_table
  attr(out, "notes") <- notes
  out
}

# Report the per-cell notes collected by WHO_mediation_analysis(): skipped
# groups and dropped indicators as messages, empty fits as one warning of
# class poemed_empty_fit (the class the single-fit warning carries, so a
# handler written for one catches the other), and failed fits as a warning
# quoting the error. Not exported.
#' @keywords internal
#' @noRd
.WHO_report_notes <- function(notes) {
  if (is.null(notes) || !nrow(notes)) return(invisible(NULL))
  label <- function(df) paste(sprintf("%s (%s)", df$group, df$detail),
                              collapse = "; ")
  count <- function(df) sprintf("%d group%s", nrow(df),
                                if (nrow(df) == 1L) "" else "s")
  skipped <- notes[notes$stage == "design", , drop = FALSE]
  missing <- notes[notes$stage == "missing", , drop = FALSE]
  empty   <- notes[notes$stage == "empty", , drop = FALSE]
  failed  <- notes[notes$stage == "fit", , drop = FALSE]
  if (nrow(skipped))
    message(sprintf("%s skipped and reported as NA: %s",
                    count(skipped), label(skipped)))
  if (nrow(missing))
    message(sprintf(paste0(
      "Indicators with missing or non-finite values were dropped before ",
      "fitting in %s: %s. Complete those values in `data` to test every ",
      "group on the full set of indicators."), count(missing), label(missing)))
  if (nrow(empty))
    warning(structure(
      class = c("poemed_empty_fit", "warning", "condition"),
      list(message = sprintf(paste0(
        "The penalized fit that HBIC selected contains no mediator in %s: ",
        "%s. The benchmark and power-enhanced tests were not computed there; ",
        "both are reported with p-values of 1, POEMED's reporting convention ",
        "for an empty selection. An empty fit is expected when no mediator ",
        "is active, but a p-value of 1 here is not evidence for the null; ",
        "see the `lambda_grid` argument of ?pe_mediation."),
        count(empty), label(empty)),
        call = NULL)))
  if (nrow(failed))
    warning(sprintf("The fit failed in %s, reported as NA: %s",
                    count(failed), label(failed)), call. = FALSE)
  invisible(NULL)
}
