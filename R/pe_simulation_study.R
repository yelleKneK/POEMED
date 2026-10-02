#' Run the Power Study Under Several Mediation Patterns at Once
#'
#' A convenience wrapper around [pe_power_curve()] that runs the size and
#' power study under more than one mediation pattern in a single call and
#' stacks the results, so a homogeneous and a contrasting curve (the two
#' panels of a figure in the article) come back together. Each pattern is
#' swept over the same signal grid.
#'
#' @details
#' Everything [pe_power_curve()] says about Monte Carlo error, empty fits,
#' and the tuning grid applies here; the article's logistic and Poisson
#' figures fix `c2` at 1 and 0.4. The patterns take their data from one continuing
#' random stream, so no two patterns share their draws. With `cores` above
#' 1 the points of `c1_grid` within a pattern do share them (see `cores` on
#' [pe_power_curve()]).
#'
#' @param n,p Observations and candidate mediators per simulated data set.
#' @param outcome Outcome type: `"continuous"`, `"binary"`, or `"count"`.
#' @param patterns Character vector of mediation patterns to run, each
#'   named once. The default runs `"homogeneous"` and `"contrasting"`.
#' @param c1_grid,c2,n_rep,alpha_level,method,error_level,lambda_grid,cores
#'   Passed to [pe_power_curve()].
#' @param seed Optional integer seed. When supplied it is set for the
#'   duration of the call and the caller's random number generator state is
#'   restored on exit. `NULL` (the default) sets no seed: the draws come from
#'   the session's random number stream, which the call advances as any
#'   random function does.
#'
#' @return A tidy `data.frame` (class `poemed_tbl`) with a leading `pattern`
#'   column and the [pe_power_curve()] columns (`c1`, `rejection_hdmm`,
#'   `rejection_pe`, `n_valid`, `n_empty`) for each pattern. The settings
#'   are recorded in the attributes `outcome`, `n`, `p`, `n_rep`, `alpha_level`,
#'   `error_level`, and `lambda_grid`.
#'
#' @inheritSection POEMED-package How to Cite
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_power_curve()] for a single pattern, [plot.poemed_tbl()] to
#'   draw one pattern's curves.
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
#' pe_simulation_study(n = 120, p = 50, outcome = "continuous",
#'                     c1_grid = c(0, 1), n_rep = 5,
#'                     lambda_grid = seq(0.05, 2, length.out = 20))
#'
#' @export
pe_simulation_study <- function(n, p,
                                outcome = c("continuous", "binary", "count"),
                                patterns = c("homogeneous", "contrasting"),
                                c1_grid = seq(0, 1, by = 0.25), c2 = 0.5,
                                n_rep = 100, alpha_level = 0.05,
                                method = c("Bonferroni", "BH", "BY"),
                                error_level = 0.05,
                                lambda_grid = seq(0.05, 10, length.out = 100),
                                cores = 1L, seed = NULL) {
  outcome <- .pe_match_choice(outcome, c("continuous", "binary", "count"),
                              "outcome")
  method <- .pe_match_choice(method, c("Bonferroni", "BH", "BY"), "method")
  patterns <- .pe_match_choice(patterns, c("homogeneous", "contrasting",
                                           "heterogeneous"), "patterns",
                               several_ok = TRUE)
  .validate_mc_design(n, p, n_rep, c1_grid, c2, cores)
  .validate_alpha_level(alpha_level, "alpha_level")
  .validate_error_level(error_level)
  # Check every pattern's design before the first one runs, so a design that
  # fails for a later pattern does not stop the study after the earlier ones.
  for (pat in patterns) .pe_mc_check_design(n, p, outcome, pat, c1_grid, c2)

  # Seed once here (rather than per pattern) so the whole study is reproducible
  # and the patterns are not handed the same stream.
  .poemed_local_seed(seed)

  parts <- lapply(patterns, function(pat) {
    pc <- pe_power_curve(n = n, p = p, outcome = outcome, pattern = pat,
                         c1_grid = c1_grid, c2 = c2, n_rep = n_rep,
                         alpha_level = alpha_level, method = method,
                         error_level = error_level, lambda_grid = lambda_grid,
                         cores = cores)
    cbind(pattern = pat, as.data.frame(unclass(pc)), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out <- .as_poemed_tbl(out)
  attr(out, "outcome")     <- outcome
  attr(out, "n")           <- n
  attr(out, "p")           <- p
  attr(out, "n_rep")       <- n_rep
  attr(out, "alpha_level") <- alpha_level
  attr(out, "error_level") <- error_level
  attr(out, "lambda_grid") <- lambda_grid
  out
}
