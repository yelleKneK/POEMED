#' Monte Carlo Size and Power Curve for the PE Mediation Tests
#'
#' Runs the article's simulation design: for each value of the
#' signal-strength scale \eqn{c_1}{c1} on a grid, it simulates many data
#' sets, applies both the benchmark Wald test and the power-enhanced test,
#' and returns their empirical rejection rates. At \eqn{c_1 = 0}{c1 = 0}
#' (the global null of no mediation) the rejection rate estimates the Type
#' I error rate; at nonzero \eqn{c_1}{c1} it estimates power. The headline
#' finding is visible directly in the table: under a contrasting pattern
#' the power-enhanced test climbs toward one as \eqn{|c_1|}{|c1|} grows,
#' while the benchmark test, whose target is a total indirect effect that
#' cancels, gains power much more slowly.
#'
#' @details
#' A rejection rate from `n_rep` replications carries a Monte Carlo
#' standard error of about \eqn{\sqrt{r(1 - r) / n_{rep}}}{sqrt(r (1 - r) /
#' n_rep)}; at 100 replications a rate near 0.05 is known to about 0.02
#' and a rate near 0.5 to about 0.05. The article uses 1000 replications. A
#' replication whose penalized fit selects no mediator counts as a
#' non-rejection (both p-values are 1), POEMED's reporting convention (the
#' article states only that \eqn{J_m = 0}{J_m = 0} when the selection is
#' empty); the number of such replications is reported in `n_empty`. A
#' replication whose fit fails outright is dropped from the denominator and
#' reported in `n_valid`, with a warning that quotes the first error. A design
#' that cannot be simulated (for example a `p` too small for the pattern, or a
#' misspelled name in `outcome_args`) stops before any replication runs, with
#' the same message at any value of `cores`.
#'
#' The tuning grid matters for what these curves show. Every replication
#' searches the same grid: the package default `seq(0.05, 10, length.out =
#' 100)`, or the grid passed as `lambda_grid`. A study of many replications
#' is where a grid tuned to the situation pays: fewer values, or a range
#' narrowed to where the HBIC minimum has been seen to fall, make each fit
#' faster at the price of a coarser choice, and whatever the grid, each
#' fit keeps the best value the grid offers. The single-fit warning for a
#' choice at an end of the grid (class `poemed_grid_boundary`) is not
#' raised inside a study, which searches a fixed grid by design. The
#' reproduction vignette states the grids of the article's own scripts.
#'
#' Under contrasting mediation the benchmark can have power although the
#' total indirect effect is zero: a penalized fit that keeps only part of
#' the canceling set leaves a total indirect effect, over the kept
#' mediators, that no longer cancels, and how often that happens depends
#' on the tuning grid. The article's logistic and Poisson figures fix
#' \eqn{c_2}{c2} at 1 and 0.4, not at this function's default of 0.5, so
#' pass `c2` to match them.
#'
#' @param n,p Number of observations and candidate mediators per data set.
#' @param outcome Outcome type passed to [simulate_mediation_data()]:
#'   `"continuous"`, `"binary"`, or `"count"`.
#' @param pattern Mediation pattern: `"homogeneous"` or `"contrasting"`
#'   (`"heterogeneous"` is a synonym for `"contrasting"`).
#' @param c1_grid Numeric vector of signal-strength scales to sweep.
#'   Default `seq(0, 1, by = 0.25)`. Include 0 to estimate the Type I
#'   error rate; the article also uses negative values.
#' @param c2 Direct effect used in the simulation. Default 0.5, the value
#'   of the article's linear simulations; its logistic figures use 1 and
#'   its Poisson figures 0.4.
#' @param n_rep Number of Monte Carlo replications per grid point. Default
#'   100. The article uses 1000.
#' @param alpha_level Significance level for the global test. A replication
#'   counts as a rejection when its p-value is at most `alpha_level`. Default
#'   0.05.
#' @param method,error_level Passed to [pe_mediation()]: the multiplicity
#'   method and target error rate for the mediator-screening step.
#' @param lambda_grid Passed to [pe_mediation()]: the tuning grid for the
#'   penalized fit, by default `seq(0.05, 10, length.out = 100)` (see
#'   Details). A shorter or narrower grid speeds a study up.
#' @param outcome_args A named list of further arguments forwarded to
#'   [simulate_mediation_data()] (for example `rho`, `q`, `d`, `tau`,
#'   `alpha_m`). Default empty. Each name must identify one argument of
#'   [simulate_mediation_data()]; an unknown or ambiguous name stops with an
#'   error rather than reaching the simulator. The design arguments this
#'   function sets itself (`n`, `p`, `outcome`, `pattern`, `c1`, `c2`) and
#'   `seed` may not appear here, whether spelled in full or abbreviated; a
#'   `seed` inside `outcome_args` would make every replication draw the same
#'   data.
#' @param cores Number of CPU cores for the Monte Carlo replications.
#'   Values above 1 fork via the base \pkg{parallel} package (Unix only).
#'   With parallel cores the replications run on the `"L'Ecuyer-CMRG"`
#'   generator (the caller's generator kind is restored on exit), so a
#'   seeded parallel run is reproducible across runs at the same `cores` but
#'   need not match a serial run. The parallel streams are handed out afresh
#'   from the same generator state at every point of `c1_grid`, so with
#'   `cores` above 1 each grid point replays the same `n_rep` random draws
#'   (common random numbers): the rates at neighboring grid points share
#'   their Monte Carlo noise, and the curve looks smoother than the
#'   separate errors of its points suggest. With `cores = 1` every grid
#'   point draws fresh data. Default 1.
#' @param seed Optional integer seed. When supplied it is set for the
#'   duration of the call and the caller's random number generator state is
#'   restored on exit. `NULL` (the default) sets no seed: the draws come from
#'   the session's random number stream, which the call advances as any
#'   random function does.
#' @param progress `TRUE` or `FALSE`; if `TRUE`, print a line per grid point
#'   as it completes. Default `FALSE`.
#'
#' @return A tidy `data.frame` of class `poemed_tbl` with one row per grid
#'   point and columns `c1`, `rejection_hdmm`, and `rejection_pe` (the
#'   empirical rejection rates of the benchmark and power-enhanced tests),
#'   `n_valid` (replications whose fit succeeded, the denominator of the
#'   rates), and `n_empty` (replications whose penalized fit selected no
#'   mediator, counted as non-rejections). The simulation settings are
#'   recorded in the attributes `outcome`, `pattern`, `n`, `p`, `n_rep`,
#'   `alpha_level`, `error_level`, and `lambda_grid`.
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
#' @seealso [pe_mediation()], [simulate_mediation_data()].
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
#' pe_power_curve(n = 200, p = 60, outcome = "continuous",
#'                pattern = "contrasting", c1_grid = c(0, 0.5, 1),
#'                n_rep = 5, lambda_grid = seq(0.05, 2, length.out = 20))
#'
#' @export
pe_power_curve <- function(n, p,
                           outcome = c("continuous", "binary", "count"),
                           pattern = c("homogeneous", "contrasting",
                                       "heterogeneous"),
                           c1_grid = seq(0, 1, by = 0.25), c2 = 0.5,
                           n_rep = 100, alpha_level = 0.05,
                           method = c("Bonferroni", "BH", "BY"),
                           error_level = 0.05,
                           lambda_grid = seq(0.05, 10, length.out = 100),
                           outcome_args = list(),
                           cores = 1L, seed = NULL, progress = FALSE) {
  outcome <- .pe_match_choice(outcome, c("continuous", "binary", "count"),
                              "outcome")
  pattern <- .pe_match_choice(pattern, c("homogeneous", "contrasting",
                                         "heterogeneous"), "pattern")
  method <- .pe_match_choice(method, c("Bonferroni", "BH", "BY"), "method")
  .validate_mc_design(n, p, n_rep, c1_grid, c2, cores)
  .validate_alpha_level(alpha_level, "alpha_level")
  .validate_error_level(error_level)
  .pe_check_flag(progress, "progress")
  outcome_args <- .validate_outcome_args(outcome_args)
  .validate_lambda_grid(lambda_grid, "lambda_grid")
  # Check once, before any replication, that the design can be simulated, so
  # a bad design reports its own message at any core count (inside a forked
  # worker it would come back as a bare "try-error" instead).
  .pe_mc_check_design(n, p, outcome, pattern, c1_grid, c2, outcome_args)

  .poemed_local_seed(seed, parallel = cores > 1L)

  reject_hdmm <- numeric(length(c1_grid))
  reject_pe   <- numeric(length(c1_grid))
  n_valid     <- integer(length(c1_grid))
  n_empty     <- integer(length(c1_grid))

  for (g in seq_along(c1_grid)) {
    c1 <- c1_grid[g]
    # One replication: simulate, fit both tests, return their two p-values
    # and whether the penalized selection was empty. An empty selection is a
    # non-rejection (both p-values 1) and is counted rather than warned about
    # here; a fit that errors returns NA with its message and is dropped from
    # the denominator rather than aborting the sweep.
    one_rep <- function(r) {
      dat <- do.call(simulate_mediation_data,
                     c(list(n = n, p = p, outcome = outcome, pattern = pattern,
                            c1 = c1, c2 = c2), outcome_args))
      .pe_mc_fit(dat, outcome = outcome, method = method,
                 error_level = error_level, lambda_grid = lambda_grid)
    }
    reps <- .pe_mc_lapply(seq_len(n_rep), one_rep, cores = cores,
                          where = sprintf("c1 = %.3g", c1))
    p_hdmm <- vapply(reps, function(x) x$p_hdmm, numeric(1))
    p_pe   <- vapply(reps, function(x) x$p_pe, numeric(1))
    empty  <- vapply(reps, function(x) x$empty, logical(1))
    errs   <- unlist(lapply(reps, function(x) x$error))
    ok <- !is.na(p_hdmm) & !is.na(p_pe)
    .pe_mc_report_failures(errs, n_rep, sprintf("c1 = %.3g", c1))
    n_valid[g]     <- sum(ok)
    n_empty[g]     <- sum(empty[ok])
    reject_hdmm[g] <- mean(p_hdmm[ok] <= alpha_level)
    reject_pe[g]   <- mean(p_pe[ok] <= alpha_level)
    if (progress)
      message(sprintf("c1 = %.3g: reject HDMM = %.3f, PE = %.3f (%d/%d valid, %d empty)",
                      c1, reject_hdmm[g], reject_pe[g], n_valid[g], n_rep,
                      n_empty[g]))
  }

  out <- data.frame(c1 = c1_grid, rejection_hdmm = reject_hdmm,
                    rejection_pe = reject_pe, n_valid = n_valid,
                    n_empty = n_empty, stringsAsFactors = FALSE)
  out <- .as_poemed_tbl(out)
  attr(out, "outcome")     <- outcome
  attr(out, "pattern")     <- pattern
  attr(out, "n")           <- n
  attr(out, "p")           <- p
  attr(out, "n_rep")       <- n_rep
  attr(out, "alpha_level") <- alpha_level
  attr(out, "error_level") <- error_level
  attr(out, "lambda_grid") <- lambda_grid
  out
}

# One Monte Carlo fit shared by the study functions: fit pe_mediation() on a
# simulated data set, muffling the classed empty-fit warning (counted instead)
# and the grid-boundary warning (a study searches a fixed grid by design), and
# turning an error into an NA record that carries the message. Returns a
# list with the two p-values, the empty flag, the fit (or NULL), and the error
# message (or NULL). Not exported.
#' @keywords internal
#' @noRd
.pe_mc_fit <- function(dat, outcome, method, error_level, lambda_grid, ...) {
  state <- new.env(parent = emptyenv()); state$empty <- FALSE
  fit <- tryCatch(
    withCallingHandlers(
      pe_mediation(dat$X, dat$Y, dat$M, Z = dat$Z, outcome = outcome,
                   method = method, error_level = error_level,
                   lambda_grid = lambda_grid, ...),
      poemed_empty_fit = function(w) {
        state$empty <- TRUE
        invokeRestart("muffleWarning")
      },
      poemed_grid_boundary = function(w) invokeRestart("muffleWarning")),
    error = function(e) e)
  if (inherits(fit, "error"))
    return(list(p_hdmm = NA_real_, p_pe = NA_real_, empty = FALSE, fit = NULL,
                error = conditionMessage(fit)))
  list(p_hdmm = fit$value[fit$term == "pval_hdmm"],
       p_pe = fit$value[fit$term == "pval_pe"], empty = state$empty, fit = fit,
       error = NULL)
}

# Report failed Monte Carlo replications: a warning naming the count and the
# first message when some failed, an error when every replication failed (a
# rate over zero valid replications is not a number). Not exported.
#' @keywords internal
#' @noRd
.pe_mc_report_failures <- function(errs, n_rep, where) {
  if (!length(errs)) return(invisible(TRUE))
  if (length(errs) >= n_rep)
    stop(sprintf("Every one of the %d replications at %s failed to fit; the first error was: %s",
                 n_rep, where, errs[1L]), call. = FALSE)
  warning(sprintf("%d of %d replications at %s failed to fit and were dropped; the first error was: %s",
                  length(errs), n_rep, where, errs[1L]), call. = FALSE)
  invisible(TRUE)
}

# Run the replications of one grid point through .pe_lapply() and stop, with
# the worker's own message, if any of them failed in a parallel worker. With
# cores above 1, parallel::mclapply() does not stop when a replication raises
# an error outside the fit's own tryCatch (in the simulation step, say): it
# returns the affected elements as "try-error" strings, or NULL when a worker
# died, with a generic warning ("all scheduled cores encountered errors in user
# code"), and the aggregation that follows would then fail with an unrelated
# message about `$` on an atomic vector. Those generic warnings are muffled,
# since the error raised here replaces them with the first worker's message.
# Not exported.
#' @keywords internal
#' @noRd
.pe_mc_lapply <- function(X, FUN, cores, where) {
  generic <- paste0("encountered errors? in user code|resulted in an error|",
                    "did not deliver (a result|results)")
  reps <- withCallingHandlers(
    .pe_lapply(X, FUN, cores = cores),
    warning = function(w) {
      if (cores > 1L && grepl(generic, conditionMessage(w)))
        invokeRestart("muffleWarning")
    })
  .pe_mc_worker_errors(reps, where)
  reps
}

#' @keywords internal
#' @noRd
.pe_mc_worker_errors <- function(reps, where) {
  bad <- vapply(reps, function(x) is.null(x) || inherits(x, "try-error"),
                logical(1L))
  if (!any(bad)) return(invisible(TRUE))
  first <- reps[[which(bad)[1L]]]
  msg <- if (is.null(first)) {
    "a parallel worker returned no result (it may have run out of memory)"
  } else if (inherits(attr(first, "condition"), "condition")) {
    conditionMessage(attr(first, "condition"))
  } else {
    trimws(as.character(first))
  }
  stop(sprintf("%d of %d replications at %s did not complete in the parallel workers; the first error was: %s",
               sum(bad), length(reps), where, msg), call. = FALSE)
}

# Validation shared by the Monte Carlo study functions. Not exported.
#' @keywords internal
#' @noRd
.validate_mc_design <- function(n, p, n_rep, c1_grid, c2, cores) {
  for (nm in c("n", "p", "n_rep")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 1 ||
        v != round(v))
      stop(sprintf("`%s` must be a single positive whole number.", nm),
           call. = FALSE)
  }
  if (!is.numeric(c1_grid) || length(c1_grid) < 1L || any(!is.finite(c1_grid)))
    stop("`c1_grid` must be a non-empty numeric vector of finite values.",
         call. = FALSE)
  if (!is.numeric(c2) || length(c2) != 1L || !is.finite(c2))
    stop("`c2` must be a single finite number.", call. = FALSE)
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
      cores < 1 || cores != round(cores))
    stop("`cores` must be a single positive whole number.", call. = FALSE)
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.validate_alpha_level <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1)
    stop(sprintf("`%s` must be a single number in (0, 1).", name), call. = FALSE)
  invisible(TRUE)
}

# `outcome_args` may carry only arguments of simulate_mediation_data() that the
# study does not set itself. The names are resolved first, the way R would
# match them when the list is spliced into the call: an abbreviation such as
# `se` would otherwise slip past a check on the full names and bind to `seed`,
# handing every replication the same data set. An unknown or ambiguous name,
# a missing name, and a name given twice all stop here, naming `outcome_args`,
# rather than surfacing later as "unused argument". Returns the list with its
# names spelled in full. Not exported.
#' @keywords internal
#' @noRd
.validate_outcome_args <- function(outcome_args) {
  if (!is.list(outcome_args))
    stop("`outcome_args` must be a list.", call. = FALSE)
  if (!length(outcome_args)) return(outcome_args)
  given <- names(outcome_args)
  if (is.null(given) || anyNA(given) || any(!nzchar(given)))
    stop("Every element of `outcome_args` must be named after an argument of ",
         "simulate_mediation_data(), for example `list(rho = 0.3)`.",
         call. = FALSE)
  formal <- names(formals(simulate_mediation_data))
  hit <- pmatch(given, formal, duplicates.ok = TRUE)
  if (anyNA(hit))
    stop("`outcome_args` has names that do not identify one argument of ",
         "simulate_mediation_data(): ",
         paste0("`", given[is.na(hit)], "`", collapse = ", "),
         ". Spell each name in full (see ?simulate_mediation_data).",
         call. = FALSE)
  full <- formal[hit]
  if (anyDuplicated(full))
    stop("`outcome_args` gives ",
         paste0("`", unique(full[duplicated(full)]), "`", collapse = ", "),
         " more than once; give each argument once.", call. = FALSE)
  reserved <- c("n", "p", "outcome", "pattern", "c1", "c2", "seed")
  bad <- full %in% reserved
  if (any(bad)) {
    lab <- ifelse(given[bad] == full[bad], paste0("`", full[bad], "`"),
                  paste0("`", full[bad], "` (given as `", given[bad], "`)"))
    stop("`outcome_args` may not contain ", paste(lab, collapse = ", "),
         "; the study sets those itself (and a `seed` there would make every ",
         "replication identical).", call. = FALSE)
  }
  names(outcome_args) <- full
  outcome_args
}

# Check, before the first replication, that simulate_mediation_data() accepts
# the design at every point of the signal grid. The check runs the simulator's
# own validation (.pe_sim_design(), which draws no random numbers), so it
# raises the simulator's own message and leaves the random stream untouched.
# Not exported.
#' @keywords internal
#' @noRd
.pe_mc_check_design <- function(n, p, outcome, pattern, c1_grid, c2,
                                outcome_args = list()) {
  for (c1 in unique(c1_grid))
    do.call(.pe_sim_design,
            c(list(n = n, p = p, outcome = outcome, pattern = pattern,
                   c1 = c1, c2 = c2), outcome_args))
  invisible(TRUE)
}

# Resolve an enumerated argument as match.arg() does, but stop with a message
# that names the argument the user passed (match.arg() calls every argument
# `arg`). An argument left at its default arrives as the whole `choices`
# vector and resolves to its first element (to all of them when
# `several_ok = TRUE`); NULL resolves to the first element, as in
# match.arg(). Unique abbreviations are accepted. With `several_ok = TRUE`
# every element must identify a choice and no choice may be named twice (a
# repeated method would be scored twice). Not exported.
#' @keywords internal
#' @noRd
.pe_match_choice <- function(arg, choices, name, several_ok = FALSE) {
  if (is.null(arg)) return(choices[1L])
  if (identical(arg, choices)) return(if (several_ok) arg else arg[1L])
  one_of <- paste0("\"", choices, "\"", collapse = ", ")
  ok <- is.character(arg) && length(arg) >= 1L && !anyNA(arg) &&
    (several_ok || length(arg) == 1L)
  hit <- if (ok) pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE) else 0L
  if (!ok || any(hit == 0L))
    stop(sprintf("`%s` must be %s %s.", name,
                 if (several_ok) "one or more of" else "one of", one_of),
         call. = FALSE)
  out <- choices[hit]
  if (anyDuplicated(out))
    stop(sprintf("`%s` names %s more than once; list each once.", name,
                 paste0("\"", unique(out[duplicated(out)]), "\"",
                        collapse = ", ")), call. = FALSE)
  out
}

# A logical flag must be a single TRUE or FALSE; NA, a string, or a vector
# would otherwise be read silently by if() or isTRUE(). Not exported.
#' @keywords internal
#' @noRd
.pe_check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  invisible(TRUE)
}
