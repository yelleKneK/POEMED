#' Monte Carlo Study of Individual-Mediator Identification (FWER / FDR)
#'
#' Evaluates how well the power-enhanced screen recovers the *individual*
#' active mediators, the experiment reported in the article's supplement.
#' For each signal-strength scale \eqn{c_1}{c1} on a grid it simulates many
#' data sets, selects a set of mediators under each multiplicity method, and
#' scores those selections against the known truth, returning the empirical
#' familywise error rate, false discovery rate, precision, and recall. This
#' complements [pe_power_curve()], which evaluates the *global* test rather
#' than which mediators are flagged.
#'
#' @details
#' The supplement of Yu and Kelley (in press) reports four metrics,
#' empirical FWER, empirical FDR, precision, and recall, without defining
#' them; this function computes them as follows. In each replication the
#' selected set is compared with the truth set. The familywise error rate
#' is the proportion of replications selecting at least one mediator
#' outside the truth set. The false discovery rate is the mean false
#' discovery proportion (zero when nothing is selected). Precision is the
#' mean proportion of selected mediators that are true (scored zero when
#' nothing is selected), and recall is the mean proportion of true
#' mediators selected. The article uses 1000 replications.
#'
#' The familywise error rate is a proportion of replications, so its Monte
#' Carlo standard error is \eqn{\sqrt{r(1 - r) / n_{valid}}}{sqrt(r (1 - r)
#' / n_valid)}, with \eqn{r}{r} the reported rate. The other three metrics
#' are means of per-replication fractions between 0 and 1, for which that
#' formula is only an upper bound. Their standard error is the standard
#' deviation of the per-replication values divided by
#' \eqn{\sqrt{n_{valid}}}{sqrt(n_valid)}, and for recall it can be as little
#' as a third of the bound.
#'
#' What counts as a true mediator is set by `truth`. Under `truth =
#' "outcome_effect"` (the default) a mediator is true when its outcome
#' coefficient \eqn{\alpha_{m,j}}{alpha_m[j]} is nonzero. This is the
#' convention of the supplement's identification tables: their
#' \eqn{c_1 = 0}{c1 = 0} rows report positive precision and recall, which
#' only a truth set that is nonempty when the exposure-on-mediator paths
#' are all zero allows. There the familywise error rate is the probability
#' of selecting a mediator with no outcome effect. Under `truth =
#' "mediation"` a mediator is true only when both paths are nonzero. That
#' is the definition of the true active set in Theorem 3 of the article,
#' against which its familywise error guarantee is stated, and of the
#' simulator's `active_mediators`. At \eqn{c_1 = 0}{c1 = 0} no mediator is
#' then active, any selection is a false positive, and recall is `NA`. For
#' nonzero \eqn{c_1}{c1} the two definitions agree under the article's
#' designs, whose exposure-on-mediator loadings are all nonzero.
#'
#' The tuning grid governs these rates as much as the screen does. Every
#' replication searches the same grid: the package default
#' `seq(0.05, 10, length.out = 100)`, or the grid passed as `lambda_grid`,
#' and each fit keeps the best value the grid offers. A shorter or narrower
#' grid speeds a study up; the single-fit warning for a choice at an end of
#' the grid is not raised inside a study, which searches a fixed grid by
#' design. The reproduction vignette states the grids of the article's own
#' scripts.
#'
#' @param n,p Number of observations and candidate mediators per data set.
#' @param outcome Outcome type passed to [simulate_mediation_data()]:
#'   `"continuous"`, `"binary"`, or `"count"`.
#' @param pattern Mediation pattern: `"homogeneous"` or `"contrasting"`
#'   (`"heterogeneous"` is a synonym for `"contrasting"`).
#' @param c1_grid Numeric vector of signal-strength scales to sweep.
#'   Default `seq(0, 1, by = 0.25)`. Include 0 to estimate the null
#'   familywise error rate.
#' @param c2 Direct effect used in the simulation. Default 0.5.
#' @param n_rep Number of Monte Carlo replications per grid point. Default
#'   200. The article uses 1000.
#' @param methods Multiplicity methods to evaluate, one or more of
#'   `"Bonferroni"` (familywise error rate), `"BH"`, and `"BY"` (false
#'   discovery rate), each named once. Default all three.
#' @param error_level Target error rate for the mediator-screening step
#'   (the FWER level for Bonferroni, the FDR level for BH and BY). Default
#'   0.05.
#' @param truth Which mediators count as truly active: `"outcome_effect"`
#'   (nonzero outcome coefficient, the convention of the supplement's
#'   identification tables) or `"mediation"` (both paths nonzero, the
#'   article's definition of the true active set in its Theorem 3). See
#'   Details.
#' @param lambda_grid Passed to [pe_mediation()]: the tuning grid for the
#'   penalized fit, by default `seq(0.05, 10, length.out = 100)` (see
#'   Details). A shorter or narrower grid speeds a study up.
#' @param outcome_args A named list of further arguments forwarded to
#'   [simulate_mediation_data()] (for example `rho`, `q`, `d`, `tau`,
#'   `alpha_m`). Default empty. Each name must identify one argument of
#'   [simulate_mediation_data()]; an unknown or ambiguous name stops with an
#'   error. The design arguments this function sets itself (`n`, `p`,
#'   `outcome`, `pattern`, `c1`, `c2`) and `seed` may not appear here,
#'   whether spelled in full or abbreviated.
#' @param cores Number of CPU cores for the replications. Values above 1
#'   fork via the base \pkg{parallel} package (Unix only) and run the
#'   replications on the `"L'Ecuyer-CMRG"` generator (the caller's generator
#'   kind is restored on exit); a seeded parallel run is reproducible across
#'   runs at the same `cores` but need not match a serial run. The parallel
#'   streams are handed out afresh from the same generator state at every
#'   point of `c1_grid`, so with `cores` above 1 each grid point replays the
#'   same `n_rep` random draws (common random numbers) and the rows at
#'   different grid points share their Monte Carlo noise. With `cores = 1`
#'   every grid point draws fresh data. Default 1.
#' @param seed Optional integer seed. When supplied it is set for the
#'   duration of the call and the caller's random number generator state is
#'   restored on exit. `NULL` (the default) sets no seed: the draws come from
#'   the session's random number stream, which the call advances as any
#'   random function does.
#' @param progress `TRUE` or `FALSE`; if `TRUE`, print a line per grid
#'   point. Default `FALSE`.
#'
#' @return A tidy `data.frame` of class `poemed_tbl` with one row per
#'   (\eqn{c_1}{c1}, method) combination and columns `c1`, `method`, `fwer`,
#'   `fdr`, `precision`, `recall` (as defined in Details), `n_valid`
#'   (replications whose fit succeeded), and `n_empty` (replications whose
#'   penalized fit selected no mediator). The settings are recorded in the
#'   attributes `outcome`, `pattern`, `n`, `p`, `n_rep`, `error_level`,
#'   `truth`, and `lambda_grid`.
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
#' @seealso [pe_power_curve()] for the global test, [pe_selection()] for the
#'   per-method selected set of a single fit.
#'
#' @family mediation simulation
#'
#' @examples
#' # n_rep = 5 keeps this example fast; a rate from five replications has a
#' # Monte Carlo standard error of up to 0.22, so read the shape, not the
#' # numbers. A reported study uses the article's design (n = 300, p = 500)
#' # and 1000 replications. A study is also where a shorter tuning grid
#' # pays: every replication searches the whole grid, so this one passes 20
#' # values from 0.05 to 2 instead of the default 100 values from 0.05 to 10
#' # (see the lambda_grid argument of ?pe_mediation).
#' set.seed(113)
#' pe_identification_study(n = 200, p = 60, outcome = "continuous",
#'                         pattern = "contrasting", c1_grid = c(0, 1),
#'                         n_rep = 5,
#'                         lambda_grid = seq(0.05, 2, length.out = 20))
#'
#' @export
pe_identification_study <- function(n, p,
                                    outcome = c("continuous", "binary",
                                                "count"),
                                    pattern = c("homogeneous", "contrasting",
                                                "heterogeneous"),
                                    c1_grid = seq(0, 1, by = 0.25), c2 = 0.5,
                                    n_rep = 200,
                                    methods = c("Bonferroni", "BH", "BY"),
                                    error_level = 0.05,
                                    truth = c("outcome_effect", "mediation"),
                                    lambda_grid = seq(0.05, 10, length.out = 100),
                                    outcome_args = list(),
                                    cores = 1L, seed = NULL, progress = FALSE) {
  outcome <- .pe_match_choice(outcome, c("continuous", "binary", "count"),
                              "outcome")
  pattern <- .pe_match_choice(pattern, c("homogeneous", "contrasting",
                                         "heterogeneous"), "pattern")
  methods <- .pe_match_choice(methods, c("Bonferroni", "BH", "BY"), "methods",
                              several_ok = TRUE)
  truth <- .pe_match_choice(truth, c("outcome_effect", "mediation"), "truth")
  .validate_mc_design(n, p, n_rep, c1_grid, c2, cores)
  .validate_error_level(error_level)
  .pe_check_flag(progress, "progress")
  outcome_args <- .validate_outcome_args(outcome_args)
  .validate_lambda_grid(lambda_grid, "lambda_grid")
  # Check the design once, before any replication, so a bad design reports
  # its own message at any core count.
  .pe_mc_check_design(n, p, outcome, pattern, c1_grid, c2, outcome_args)

  .poemed_local_seed(seed, parallel = cores > 1L)

  # The selected set under each method, for one fit. report_all_methods
  # populates selection_by_method; the empty-fit path leaves it NULL (no
  # mediator selected under any method).
  method_sets <- function(fit) {
    sbm <- attr(fit, "selection_by_method")
    if (is.null(sbm)) {
      act <- attr(fit, "selected_mediators")
      sbm <- stats::setNames(rep(list(integer(0)), length(methods)), methods)
      prim <- attr(fit, "method")
      if (prim %in% methods) sbm[[prim]] <- act
    }
    sbm[methods]
  }
  na_row <- function() data.frame(method = methods, any_fp = NA_real_,
                                  fdp = NA_real_, precision = NA_real_,
                                  recall = NA_real_, valid = FALSE,
                                  empty = FALSE, stringsAsFactors = FALSE)

  res <- vector("list", length(c1_grid))
  for (g in seq_along(c1_grid)) {
    c1 <- c1_grid[g]
    # One replication: simulate, fit with all requested methods, and score each
    # method's selection against the truth set. A failed fit returns a row of
    # NAs flagged invalid and is dropped from the denominators; an empty
    # penalized selection is scored (nothing selected) and counted.
    one_rep <- function(r) {
      dat <- do.call(simulate_mediation_data,
                     c(list(n = n, p = p, outcome = outcome, pattern = pattern,
                            c1 = c1, c2 = c2), outcome_args))
      true_set <- if (truth == "outcome_effect") which(dat$alpha_m != 0)
                  else dat$active_mediators
      mc <- .pe_mc_fit(dat, outcome = outcome, method = methods[1L],
                       error_level = error_level, lambda_grid = lambda_grid,
                       report_all_methods = length(methods) > 1L)
      if (is.null(mc$fit)) {
        out <- na_row(); attr(out, "error") <- mc$error; return(out)
      }
      sets <- method_sets(mc$fit)
      do.call(rbind, lapply(methods, function(m) {
        sel <- sets[[m]]
        tp <- length(intersect(sel, true_set))
        fp <- length(setdiff(sel, true_set))
        n_sel <- length(sel); n_tru <- length(true_set)
        data.frame(
          method = m,
          any_fp = as.numeric(fp >= 1L),
          fdp = if (n_sel > 0L) fp / n_sel else 0,
          precision = if (n_sel > 0L) tp / n_sel else 0,
          recall = if (n_tru > 0L) tp / n_tru else NA_real_,
          valid = TRUE, empty = mc$empty, stringsAsFactors = FALSE)
      }))
    }
    rep_list <- .pe_mc_lapply(seq_len(n_rep), one_rep, cores = cores,
                              where = sprintf("c1 = %.3g", c1))
    errs <- unlist(lapply(rep_list, function(x) attr(x, "error")))
    .pe_mc_report_failures(errs, n_rep, sprintf("c1 = %.3g", c1))
    reps <- do.call(rbind, rep_list)

    # Aggregate over the valid replications, per method.
    avg <- function(x) if (any(!is.na(x))) mean(x, na.rm = TRUE) else NA_real_
    rows <- lapply(methods, function(m) {
      sub <- reps[reps$method == m & reps$valid, , drop = FALSE]
      data.frame(c1 = c1, method = m,
                 fwer = avg(sub$any_fp), fdr = avg(sub$fdp),
                 precision = avg(sub$precision), recall = avg(sub$recall),
                 n_valid = nrow(sub), n_empty = sum(sub$empty),
                 stringsAsFactors = FALSE)
    })
    res[[g]] <- do.call(rbind, rows)
    if (progress)
      message(sprintf("c1 = %.3g done (%d methods x %d reps)",
                      c1, length(methods), n_rep))
  }

  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out <- .as_poemed_tbl(out)
  attr(out, "outcome")     <- outcome
  attr(out, "pattern")     <- pattern
  attr(out, "n")           <- n
  attr(out, "p")           <- p
  attr(out, "n_rep")       <- n_rep
  attr(out, "error_level") <- error_level
  attr(out, "truth")       <- truth
  attr(out, "lambda_grid") <- lambda_grid
  out
}
