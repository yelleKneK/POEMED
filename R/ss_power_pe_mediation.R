#' Simulation-Based Sample Size Planning for the Power-Enhanced Mediation Test
#'
#' Estimates, by Monte Carlo simulation, the sample size needed for the
#' power-enhanced test to reach a target power under a given mediation
#' pattern and signal strength. It evaluates the empirical power of the PE
#' test (and, for comparison, the benchmark Wald test) at each candidate
#' sample size on a grid, and reports the smallest grid value that reaches
#' the target. This is the design counterpart of the analysis function
#' [pe_mediation()]: use it to plan a study, then [pe_mediation()] to
#' analyze it.
#'
#' @details
#' This is a simulation-based procedure, not a closed-form power formula:
#' no analytic power expression exists for the PE test, so sample size is
#' instead determined by direct Monte Carlo simulation. Power is estimated
#' by simulating `n_rep` data sets at each candidate
#' `n` with [simulate_mediation_data()] and recording how often the test
#' rejects at level `alpha_level`. Because the estimate is a Monte Carlo
#' proportion, it carries simulation error of roughly
#' \eqn{\sqrt{power(1 - power) / n\_rep}}{sqrt(power (1 - power) / n_rep)};
#' raise `n_rep` for a smoother curve and a more stable recommendation. The
#' grid approach (rather than a root search) is deliberate: the power curve
#' is monotone but noisy, so a search can stop early on a lucky draw,
#' whereas the grid shows the whole trajectory.
#'
#' An estimated power from `n_rep` replications moves in steps of
#' `1 / n_rep`, so a target above `1 - 1 / n_rep` can be reached only when
#' every replication rejects, which is weak evidence that the true power
#' reaches it. Such a target draws a warning of class
#' `poemed_target_unresolvable` naming the smallest `n_rep` that resolves
#' it. When no candidate sample size reaches the target, the recommendation
#' is `NA` and a warning of class `poemed_target_not_reached` gives the
#' largest estimated power and where it occurred, so the grid can be
#' extended. The printed table ends with the target and the recommendation
#' (or the statement that the grid did not reach the target).
#'
#' @param outcome Outcome type: `"continuous"`, `"binary"`, or `"count"`.
#' @param pattern Mediation pattern: `"homogeneous"` or `"contrasting"`
#'   (`"heterogeneous"` is a synonym for `"contrasting"`).
#' @param p Number of candidate mediators.
#' @param n_grid Vector of candidate sample sizes to evaluate: whole numbers
#'   of at least 10 (a value listed twice is evaluated once). They are
#'   evaluated, and reported, in
#'   increasing order.
#' @param c1 Signal strength (the exposure-on-mediator scale), a single
#'   finite nonzero number; at `c1 = 0` the rejection rate is the test's
#'   size, which no sample size raises to a power target. Default 1.
#' @param c2 Direct effect. Default 0.5, the value of the article's linear
#'   simulations; its logistic figures use 1 and its Poisson figures 0.4.
#' @param target_power Desired power for the PE test, a single number
#'   between 0 and 1 (see Details for the replications it needs). Default
#'   0.8.
#' @param n_rep Monte Carlo replications per candidate `n`. Default 200.
#' @param alpha_level Significance level. Default 0.05.
#' @param method,error_level Passed to [pe_mediation()].
#' @param lambda_grid Passed to [pe_mediation()]: the tuning grid searched
#'   at each candidate `n`, by default `seq(0.05, 10, length.out = 100)`. A
#'   plan runs `n_rep` fits at every candidate `n`, so a shorter or narrower
#'   grid speeds it up.
#' @param cores Number of CPU cores (Unix forking). A seeded parallel run is
#'   reproducible across runs at the same `cores` but need not match a
#'   serial run. Default 1.
#' @param seed Optional integer seed. When supplied it is set for the
#'   duration of the call and the caller's random number generator state is
#'   restored on exit. `NULL` (the default) sets no seed: the draws come from
#'   the session's random number stream, which the call advances as any
#'   random function does.
#'
#' @return A tidy `data.frame` (class `poemed_ss_plan`, a `poemed_tbl`) with
#'   one row per candidate sample size and columns `n`, `power_pe`,
#'   `power_hdmm`, `n_valid`, and `n_empty` (see [pe_power_curve()]). The
#'   smallest `n` reaching `target_power` for the PE test is stored in the
#'   `"recommended_n"` attribute (`NA` if no grid value reaches it), with
#'   `"target_power"`, `"outcome"`, `"pattern"`, `"c1"`, `"c2"`, `"n_rep"`,
#'   and `"alpha_level"`. Printing the table adds a line with the target and
#'   the recommendation.
#'
#' @inheritSection POEMED-package How to Cite
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_power_curve()], [pe_mediation()],
#'   [simulate_mediation_data()].
#'
#' @family mediation simulation
#'
#' @examples
#' # n_rep = 5 keeps this example fast; a power estimate from five
#' # replications has a Monte Carlo standard error of up to 0.22, so the
#' # recommended n here is a demonstration, not a plan. Planning a study
#' # deserves n_rep of several hundred. A plan is also where a shorter
#' # tuning grid pays: every replication searches the whole grid, so this
#' # one passes 20 values from 0.05 to 2 instead of the default 100 values
#' # from 0.05 to 10 (see the lambda_grid argument of ?pe_mediation).
#' set.seed(113)
#' plan <- ss_power_pe_mediation(outcome = "continuous", pattern = "contrasting",
#'                               p = 50, n_grid = c(150, 250), n_rep = 5,
#'                               lambda_grid = seq(0.05, 2, length.out = 20))
#' plan                        # the last line gives the recommendation
#' attr(plan, "recommended_n") # 150, the smallest grid n reaching 0.8
#'
#' @export
ss_power_pe_mediation <- function(outcome = c("continuous", "binary", "count"),
                                  pattern = c("homogeneous", "contrasting",
                                              "heterogeneous"),
                                  p, n_grid, c1 = 1, c2 = 0.5,
                                  target_power = 0.8, n_rep = 200, alpha_level = 0.05,
                                  method = c("Bonferroni", "BH", "BY"),
                                  error_level = 0.05,
                                  lambda_grid = seq(0.05, 10, length.out = 100),
                                  cores = 1L, seed = NULL) {
  outcome <- .pe_match_choice(outcome, c("continuous", "binary", "count"),
                              "outcome")
  pattern <- .pe_match_choice(pattern, c("homogeneous", "contrasting",
                                         "heterogeneous"), "pattern")
  method <- .pe_match_choice(method, c("Bonferroni", "BH", "BY"), "method")
  if (!is.numeric(n_grid) || length(n_grid) < 1L || any(!is.finite(n_grid)) ||
      any(n_grid < 10) || any(n_grid != round(n_grid)))
    stop("`n_grid` must be a vector of whole-number sample sizes, each at ",
         "least 10.", call. = FALSE)
  if (!is.numeric(c1) || length(c1) != 1L || !is.finite(c1))
    stop("`c1` must be a single finite number (one signal strength to plan for).",
         call. = FALSE)
  if (c1 == 0)
    stop("`c1 = 0` is the global null, where the rejection rate of the ",
         "power-enhanced test is its size rather than its power, so no ",
         "sample size reaches a power target. Pass the nonzero signal ",
         "strength to plan for.", call. = FALSE)
  .validate_mc_design(n_grid[1L], p, n_rep, c1, c2, cores)
  .validate_alpha_level(alpha_level, "alpha_level")
  .validate_error_level(error_level)
  if (!is.numeric(target_power) || length(target_power) != 1L ||
      !is.finite(target_power) || target_power <= 0 || target_power >= 1)
    stop("`target_power` must be a single number in (0, 1).", call. = FALSE)
  .validate_lambda_grid(lambda_grid, "lambda_grid")
  .pe_mc_check_design(min(n_grid), p, outcome, pattern, c1, c2)

  # A target above 1 - 1/n_rep is reached only when every replication
  # rejects (the estimate moves in steps of 1/n_rep); the tolerance keeps a
  # target that equals a step, such as 0.8 with n_rep = 5, resolvable.
  if (target_power > 1 - 1 / n_rep + 1e-9) {
    need <- ceiling(1 / (1 - target_power) - 1e-9)
    .pe_warn(sprintf(paste0(
      "`target_power = %s` cannot be resolved with `n_rep = %d`: the ",
      "estimated power moves in steps of 1/%d, so only a run in which ",
      "every replication rejects reaches the target. Increase `n_rep` to ",
      "at least %d."), format(target_power), as.integer(n_rep),
      as.integer(n_rep), as.integer(need)), "poemed_target_unresolvable")
  }

  .poemed_local_seed(seed)

  # Evaluate power at each candidate n by reusing the single-point power curve
  # (c1 fixed); the c1 = 0 size case is not needed here, so the grid is one c1.
  # A sample size listed twice is evaluated once.
  n_grid <- sort(unique(as.integer(n_grid)))
  rows <- lapply(n_grid, function(nn) {
    pc <- pe_power_curve(n = nn, p = p, outcome = outcome, pattern = pattern,
                         c1_grid = c1, c2 = c2, n_rep = n_rep, alpha_level = alpha_level,
                         method = method, error_level = error_level,
                         lambda_grid = lambda_grid, cores = cores)
    data.frame(n = nn, power_pe = pc$rejection_pe[1],
               power_hdmm = pc$rejection_hdmm[1], n_valid = pc$n_valid[1],
               n_empty = pc$n_empty[1], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)

  # Recommended n: the smallest grid value whose PE power reaches the target.
  reached <- which(out$power_pe >= target_power)
  recommended_n <- if (length(reached)) out$n[min(reached)] else NA_integer_

  out <- .as_poemed_tbl(out)
  class(out) <- c("poemed_ss_plan", class(out))
  attr(out, "recommended_n") <- recommended_n
  attr(out, "target_power")  <- target_power
  attr(out, "outcome")       <- outcome
  attr(out, "pattern")       <- pattern
  attr(out, "c1")            <- c1
  attr(out, "c2")            <- c2
  attr(out, "n_rep")         <- n_rep
  attr(out, "alpha_level")   <- alpha_level

  if (is.na(recommended_n))
    .pe_warn(paste0(.pe_ss_verdict(out),
                    " Extend `n_grid` upward to find the sample size."),
             "poemed_target_not_reached")
  out
}

# Signal a warning carrying a condition class, so a caller can catch or
# muffle it by class rather than by matching its text. Not exported.
#' @keywords internal
#' @noRd
.pe_warn <- function(message, class) {
  warning(structure(list(message = message, call = NULL),
                    class = c(class, "warning", "condition")))
}

# The one-line verdict a plan carries: the target and the recommended n, or,
# when no candidate reached the target, the largest estimated power and the n
# where it occurred. Used by the print method and the not-reached warning.
# Display only; nothing here feeds a returned value. Not exported.
#' @keywords internal
#' @noRd
.pe_ss_verdict <- function(x) {
  target <- attr(x, "target_power")
  rec <- attr(x, "recommended_n")
  n_rep <- attr(x, "n_rep")
  head <- sprintf("Target power %s for the power-enhanced test", format(target))
  if (!is.na(rec))
    return(sprintf("%s: recommended n = %s, the smallest n in the grid that reaches it (%s replications per n).",
                   head, format(rec), format(n_rep)))
  pw <- x$power_pe
  if (!length(pw) || all(is.na(pw)))
    return(sprintf("%s: not reached in the grid.", head))
  best <- which.max(pw)
  sprintf("%s: not reached in the grid; the largest estimated power is %s, at n = %s (%s replications per n).",
          head, format(pw[best], digits = 3), format(x$n[best]), format(n_rep))
}

#' @exportS3Method print poemed_ss_plan
print.poemed_ss_plan <- function(x, ...) {
  NextMethod()
  if (!is.null(attr(x, "target_power")) && !is.null(attr(x, "recommended_n")))
    cat(.pe_ss_verdict(x), "\n", sep = "")
  invisible(x)
}
