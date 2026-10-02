#' Per-Mediator Detail Behind a Power-Enhanced Mediation Test
#'
#' Returns the per-mediator table that explains a [pe_mediation()] result.
#' The table has one row for each candidate mediator, that is, each mediator
#' the penalized fit retained with a nonzero coefficient. Each row gives the
#' mediator's two path statistics, its screening p-value, and whether the
#' screen selected it. This is the "why" behind the selected set,
#' so a user can see which candidates were close to the screening threshold
#' and which were screened out.
#'
#' Two steps narrow the mediators, and the table keeps them apart. The
#' penalized fit keeps the candidates, which are the rows of the table. The
#' multiplicity screen then selects mediators among those candidates, which
#' are the rows whose `selected` column is `TRUE`. A row with
#' `selected = FALSE` is therefore a candidate that the screen did not
#' select, not a mediator the penalized fit dropped.
#'
#' @param fit A result from [pe_mediation()] or one of its workers.
#'
#' @return A tidy `data.frame` (class `poemed_tbl`) with columns:
#' \describe{
#'   \item{mediator}{Column position of the mediator in the original `M`.}
#'   \item{name}{The mediator's column name in `M`; present only when `M`
#'     has column names.}
#'   \item{t_outcome}{The standardized mediator-on-outcome statistic
#'     (the \eqn{M \to Y}{M -> Y} path, \eqn{\hat\alpha_{m,j} /
#'     \hat\sigma_{m,j}}{hat(alpha)_m[j] / hat(sigma)_m[j]}).}
#'   \item{t_exposure}{The standardized exposure-on-mediator statistic
#'     (the \eqn{X \to M}{X -> M} path). With one exposure it is the signed
#'     \emph{t} statistic, so a negative value marks a negative path. With
#'     several exposures it is the largest absolute \emph{t} statistic across
#'     the exposures, so it is never negative and does not show the path's
#'     direction; regress the mediator on the exposures and any
#'     confounders (for example with `lm()`) to see the signs.}
#'   \item{screen_p}{The screening p-value, the larger of the two path
#'     p-values (the smaller of these across exposures when there are
#'     several). A mediator is selected when this clears the multiplicity
#'     threshold.}
#'   \item{selected}{Logical; `TRUE` when the screen selected the mediator
#'     under the fit's `method`, so the `TRUE` rows are the fit's selected
#'     set. Every row is a candidate the penalized fit retained, so `FALSE`
#'     marks a candidate the screen did not select.}
#' }
#'   The table is empty when the penalized fit retained no candidate
#'   mediators.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()], [pe_selection()].
#'
#' @family mediation tests
#'
#' @examples
#' set.seed(113)
#' d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
#'                              pattern = "contrasting", c1 = 1)
#' fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
#' # HBIC chose the smallest value of the default grid here, so the fit
#' # warned; see the lambda_grid argument of ?pe_mediation.
#' pe_mediators(fit)
#'
#' @export
pe_mediators <- function(fit) {
  tab <- attr(fit, "mediator_table")
  if (is.null(tab))
    stop("`fit` does not carry a mediator table; was it produced by ",
         "pe_mediation()?", call. = FALSE)
  .as_poemed_tbl(tab)
}

#' Selected Mediators Under Each Multiplicity Method
#'
#' Summarizes which mediators a [pe_mediation()] fit selected.
#' When the fit was produced with `report_all_methods = TRUE`, this
#' reports the selected set under each of the three multiplicity methods
#' (Bonferroni for familywise error rate control, Benjamini-Hochberg and
#' Benjamini-Yekutieli for false discovery rate control), so their
#' conservativeness can be compared on the same fit. Otherwise it reports
#' the single method that was used.
#'
#' @param fit A result from [pe_mediation()] or one of its workers.
#'
#' @return A tidy `data.frame` of class `poemed_tbl` (so its p-values print
#'   to fixed decimals) with one row per multiplicity method and columns
#'   `method`, `n_selected` (number of selected mediators), `j_pe` (the power
#'   enhancement component \eqn{J_m}{J_m} under that method's screen), `stat_pe`
#'   and `pval_pe` (the resulting power-enhanced statistic \eqn{M_{PE}}{M_PE}
#'   and its p-value), and `selected_mediators` (their column names when `M` has
#'   them, otherwise their column positions, comma-separated, or `"none"`).
#'   Because each multiplicity method screens a different selected set into
#'   \eqn{J_m}{J_m}, the `pval_pe` column gives the PE / PE_BH / PE_BY global
#'   p-values side by side, as the article's extended data-analysis tables
#'   report them.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()] (its `report_all_methods` argument),
#'   [pe_mediators()].
#'
#' @family mediation tests
#'
#' @examples
#' set.seed(113)
#' d <- simulate_mediation_data(n = 300, p = 100, outcome = "continuous",
#'                              pattern = "homogeneous", c1 = 1)
#' fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
#'                     report_all_methods = TRUE)
#' # On this homogeneous design (active mediators 1 to 5) the screens
#' # differ: Bonferroni, which controls the familywise error rate, selects
#' # three mediators; BH, the false discovery rate screen for independent or
#' # positively dependent mediators, selects five, one of them inactive; BY,
#' # valid under arbitrary dependence, selects the same three as Bonferroni.
#' pe_selection(fit)
#' d$active_mediators
#' # HBIC chose 0.05, the smallest value of the default grid, so the fit
#' # warned. On a grid that extends below it and is spaced more finely,
#' # seq(0.01, 10, length.out = 200), it chooses about 0.11, and the
#' # screens then select 3 and 4 (Bonferroni) and 2, 3, and 4 (BH and BY):
#' # the selected set depends on the grid as well as on the screen.
#'
#' @export
pe_selection <- function(fit) {
  # Every fit carries its selected set and the method that produced it; an
  # object without them (a WHO table, a data frame, NULL) is not a fit.
  if (is.null(attr(fit, "selected_mediators")) || is.null(attr(fit, "method")))
    stop("`fit` does not carry a selection record; pass a result from ",
         "pe_mediation() or one of its workers.", call. = FALSE)
  by_method <- attr(fit, "selection_by_method")
  # When only one method was run, build a single-row summary from the stored
  # selected set and method label.
  if (is.null(by_method)) {
    sel <- attr(fit, "selected_mediators")
    by_method <- stats::setNames(list(sel), attr(fit, "method"))
  }
  meths <- names(by_method)
  # The per-method PE statistics travel as the pe_by_method attribute; align
  # them to the methods in by_method (NA if, for an old fit, they are absent).
  pbm <- attr(fit, "pe_by_method")
  pick <- function(col) if (is.null(pbm)) rep(NA_real_, length(meths)) else
    pbm[[col]][match(meths, pbm$method)]
  # Label the selected set by column name when `M` had names.
  nm <- attr(fit, "mediator_names")
  lab <- function(a) if (!length(a)) "none" else
    paste(if (is.null(nm)) a else nm[a], collapse = ", ")
  data.frame(
    method = meths,
    n_selected = vapply(by_method, length, integer(1)),
    j_pe = pick("j_pe"),
    stat_pe = pick("stat_pe"),
    pval_pe = pick("pval_pe"),
    selected_mediators = vapply(by_method, lab, character(1)),
    row.names = NULL, stringsAsFactors = FALSE) |> .as_poemed_tbl()
}
