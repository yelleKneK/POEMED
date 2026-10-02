# Input guards of the Monte Carlo study functions (pe_power_curve(),
# pe_identification_study(), pe_simulation_study(), ss_power_pe_mediation()).
# Every guard is made to fire here, and each message must name the argument
# the user passed. The calls stop before any replication runs, so they cost
# almost nothing.

big <- c(5, 10)   # a tuning grid that selects no mediator at these designs

# Muffle the planner's classed warnings where a test is about something else.
quiet_plan <- function(expr) {
  withCallingHandlers(expr,
    poemed_target_unresolvable = function(w) invokeRestart("muffleWarning"),
    poemed_target_not_reached = function(w) invokeRestart("muffleWarning"))
}

test_that("an abbreviated name in outcome_args cannot reach `seed`", {
  # `se` and `see` would bind to `seed` by partial matching and hand every
  # replication the same data set.
  for (nm in c("se", "see")) {
    args <- stats::setNames(list(1), nm)
    expect_error(pe_power_curve(n = 40, p = 10, c1_grid = 1, n_rep = 2,
                                outcome_args = args),
                 "`outcome_args` may not contain `seed` \\(given as")
    expect_error(pe_identification_study(n = 40, p = 10, c1_grid = 1,
                                          n_rep = 2, outcome_args = args),
                 "`outcome_args` may not contain `seed`")
  }
  # The same for the other reserved names, spelled in full or abbreviated.
  expect_error(pe_power_curve(n = 40, p = 10, n_rep = 2,
                              outcome_args = list(pa = "homogeneous")),
               "may not contain `pattern`")
  expect_error(pe_power_curve(n = 40, p = 10, n_rep = 2,
                              outcome_args = list(outcome = "binary")),
               "may not contain `outcome`")
})

test_that("outcome_args must be a named list of simulator arguments", {
  expect_error(pe_power_curve(n = 40, p = 10, n_rep = 2, outcome_args = 1),
               "`outcome_args` must be a list")
  expect_error(pe_power_curve(n = 40, p = 10, n_rep = 2,
                              outcome_args = list(2)),
               "Every element of `outcome_args` must be named")
  expect_error(pe_power_curve(n = 40, p = 10, n_rep = 2,
                              outcome_args = list(rho = 0.3, 2)),
               "Every element of `outcome_args` must be named")
  # Unknown and ambiguous names stop up front, naming outcome_args, rather
  # than as R's "unused argument" from inside a replication.
  expect_error(pe_power_curve(n = 40, p = 10, n_rep = 2,
                              outcome_args = list(foo = 1)),
               "`outcome_args` has names .*`foo`")
  expect_error(pe_identification_study(n = 40, p = 10, n_rep = 2,
                                       outcome_args = list(alpha = 1)),
               "`outcome_args` has names .*`alpha`")
  expect_error(pe_power_curve(n = 40, p = 10, n_rep = 2,
                              outcome_args = list(rho = 0.3, rh = 0.4)),
               "`outcome_args` gives `rho` more than once")
})

test_that("an unambiguous abbreviation in outcome_args is read as the full name", {
  a <- pe_power_curve(n = 60, p = 15, pattern = "contrasting", c1_grid = 1,
                      n_rep = 1, outcome_args = list(rh = 0.3), seed = 113)
  b <- pe_power_curve(n = 60, p = 15, pattern = "contrasting", c1_grid = 1,
                      n_rep = 1, outcome_args = list(rho = 0.3), seed = 113)
  expect_identical(a, b)
})

test_that("a design the simulator refuses stops before any replication, at any core count", {
  for (cores in c(1L, 2L)) {
    expect_error(pe_power_curve(n = 40, p = 3, pattern = "contrasting",
                                c1_grid = 1, n_rep = 2, cores = cores),
                 "`p` is too small for this pattern")
    expect_error(pe_power_curve(n = 40, p = 10, c1_grid = 1, n_rep = 2,
                                cores = cores,
                                outcome_args = list(rho = 2)),
                 "`rho` must lie strictly between -1 and 1")
    expect_error(pe_identification_study(n = 40, p = 10, c1_grid = 1,
                                         n_rep = 2, cores = cores,
                                         outcome_args = list(alpha_z = 1)),
                 "`alpha_z` describes confounders")
  }
  expect_error(pe_power_curve(n = 1, p = 10, n_rep = 2),
               "`n` must be a single whole number of at least 2")
  # pe_simulation_study() checks each pattern's design before running any.
  expect_error(pe_simulation_study(n = 40, p = 4, n_rep = 2, c1_grid = 1,
                                   patterns = c("contrasting", "homogeneous")),
               "`p` must be at least 5 for the default `tau`")
  expect_error(quiet_plan(ss_power_pe_mediation(p = 3, pattern = "contrasting",
                                                n_grid = 40, n_rep = 2)),
               "`p` is too small for this pattern")
})

test_that("an error inside a parallel worker reaches the user with its own message", {
  # The helper that reads mclapply()'s output: a "try-error" element carries
  # the worker's condition, and a NULL element is a worker that died.
  te <- structure("Error in f() : boom\n", class = "try-error",
                  condition = simpleError("boom"))
  expect_error(POEMED:::.pe_mc_worker_errors(list(list(p_hdmm = 1), te),
                                             "c1 = 1"),
               "1 of 2 replications at c1 = 1 did not complete .*: boom")
  expect_error(POEMED:::.pe_mc_worker_errors(list(NULL), "c1 = 1"),
               "returned no result")
  expect_error(POEMED:::.pe_mc_worker_errors(
    list(structure("Error : bare\n", class = "try-error")), "c1 = 1"),
    "Error : bare")
  expect_invisible(POEMED:::.pe_mc_worker_errors(list(list(p_hdmm = 1)),
                                                 "c1 = 1"))
  # End to end: a simulation step that fails inside forked workers.
  skip_on_cran()
  skip_on_os("windows")
  local_mocked_bindings(simulate_mediation_data = function(...)
    stop("simulated failure"), .package = "POEMED")
  expect_error(pe_power_curve(n = 40, p = 10, c1_grid = 1, n_rep = 2,
                              cores = 2),
               "did not complete in the parallel workers; the first error was: simulated failure")
  expect_error(pe_identification_study(n = 40, p = 10, c1_grid = 1,
                                       n_rep = 2, cores = 2),
               "the first error was: simulated failure")
})

test_that("`progress` must be TRUE or FALSE", {
  for (v in list(NA, "yes", c(TRUE, FALSE))) {
    expect_error(pe_power_curve(n = 40, p = 10, n_rep = 2, c1_grid = 0,
                                progress = v),
                 "`progress` must be TRUE or FALSE")
    expect_error(pe_identification_study(n = 40, p = 10, n_rep = 2,
                                         c1_grid = 0, progress = v),
                 "`progress` must be TRUE or FALSE")
  }
})

test_that("enumerated arguments name themselves when they do not match", {
  expect_error(pe_power_curve(n = 40, p = 10, n_rep = 2, outcome = "xyz"),
               "`outcome` must be one of \"continuous\", \"binary\", \"count\"")
  expect_error(pe_power_curve(n = 40, p = 10, n_rep = 2, pattern = NA),
               "`pattern` must be one of")
  expect_error(pe_power_curve(n = 40, p = 10, n_rep = 2, method = "holm"),
               "`method` must be one of \"Bonferroni\", \"BH\", \"BY\"")
  expect_error(pe_power_curve(n = 40, p = 10, n_rep = 2,
                              pattern = c("homogeneous", "contrasting")),
               "`pattern` must be one of")
  expect_error(pe_identification_study(n = 40, p = 10, n_rep = 2,
                                       methods = c("BH", "holm")),
               "`methods` must be one or more of")
  expect_error(pe_identification_study(n = 40, p = 10, n_rep = 2,
                                       truth = "both"),
               "`truth` must be one of \"outcome_effect\", \"mediation\"")
  expect_error(pe_simulation_study(n = 40, p = 10, n_rep = 2,
                                   patterns = "flat"),
               "`patterns` must be one or more of")
  expect_error(pe_simulation_study(n = 40, p = 10, n_rep = 2, outcome = 1),
               "`outcome` must be one of")
  expect_error(ss_power_pe_mediation(p = 10, n_grid = 40, n_rep = 2,
                                     method = "BHY"),
               "`method` must be one of")
  expect_error(ss_power_pe_mediation(outcome = "gaussian", p = 10,
                                     n_grid = 40, n_rep = 2),
               "`outcome` must be one of")
})

test_that("a method or pattern listed twice is refused, not scored twice", {
  expect_error(pe_identification_study(n = 40, p = 10, c1_grid = 1, n_rep = 2,
                                       methods = c("BH", "BH")),
               "`methods` names \"BH\" more than once")
  expect_error(pe_identification_study(n = 40, p = 10, c1_grid = 1, n_rep = 2,
                                       methods = c("BH", "B")),
               "`methods` must be one or more of")
  expect_error(pe_simulation_study(n = 40, p = 10, n_rep = 2,
                                   patterns = c("contrasting", "contrasting")),
               "`patterns` names \"contrasting\" more than once")
})

test_that("the planner validates every element of n_grid", {
  expect_error(ss_power_pe_mediation(p = 10, n_grid = c(40, 60.7), n_rep = 2),
               "`n_grid` must be a vector of whole-number sample sizes")
  expect_error(ss_power_pe_mediation(p = 10, n_grid = c(40.5, 60), n_rep = 2),
               "`n_grid` must be")
  expect_error(ss_power_pe_mediation(p = 10, n_grid = c(5, 60), n_rep = 2),
               "each at least 10")
  expect_error(ss_power_pe_mediation(p = 10, n_grid = c(40, NA), n_rep = 2),
               "`n_grid` must be")
})

test_that("a sample size listed twice in n_grid is evaluated once, quietly", {
  plan <- function(n_grid)
    quiet_plan(ss_power_pe_mediation(p = 10, n_grid = n_grid, n_rep = 2,
                                     lambda_grid = big, seed = 113))
  expect_no_message(dup <- plan(c(60, 40, 60)))
  expect_identical(dup$n, c(40L, 60L))
  expect_identical(dup, plan(c(40, 60)))
})

test_that("the planner checks its tuning grid before any replication runs", {
  expect_error(ss_power_pe_mediation(p = 10, n_grid = 40, n_rep = 2,
                                     lambda_grid = -1),
               "`lambda_grid` must be a non-empty numeric vector with all elements > 0.",
               fixed = TRUE)
  expect_error(ss_power_pe_mediation(p = 10, n_grid = 40, n_rep = 2,
                                     lambda_grid = numeric(0)),
               "`lambda_grid` must be", fixed = TRUE)
  expect_error(ss_power_pe_mediation(p = 10, n_grid = 40, n_rep = 2,
                                     lambda_grid = "a"),
               "`lambda_grid` must be a non-empty numeric vector", fixed = TRUE)
})

test_that("the planner refuses c1 = 0, where power is size", {
  expect_error(ss_power_pe_mediation(p = 10, n_grid = 40, n_rep = 2, c1 = 0),
               "`c1 = 0` is the global null")
})

test_that("a target the replication count cannot resolve draws a classed warning", {
  # With n_rep = 5 the estimate moves in steps of 0.2, so 0.8 is resolvable
  # and 0.9 is not; 0.9 needs n_rep of at least 10.
  expect_warning(
    quiet_target <- withCallingHandlers(
      ss_power_pe_mediation(p = 10, n_grid = 40, n_rep = 5, target_power = 0.9,
                            lambda_grid = big, seed = 113),
      poemed_target_not_reached = function(w) invokeRestart("muffleWarning")),
    class = "poemed_target_unresolvable")
  w <- tryCatch(
    ss_power_pe_mediation(p = 10, n_grid = 40, n_rep = 5, target_power = 0.9,
                          lambda_grid = big, seed = 113),
    poemed_target_unresolvable = function(w) w)
  expect_match(conditionMessage(w), "`target_power = 0.9` cannot be resolved with `n_rep = 5`")
  expect_match(conditionMessage(w), "Increase `n_rep` to at least 10")
  # A target equal to a step is resolvable: no such warning.
  seen <- new.env(parent = emptyenv()); seen$classes <- character()
  withCallingHandlers(
    ss_power_pe_mediation(p = 10, n_grid = 40, n_rep = 5, target_power = 0.8,
                          lambda_grid = big, seed = 113),
    warning = function(w) {
      seen$classes <- c(seen$classes, class(w)[1L])
      invokeRestart("muffleWarning")
    })
  expect_false("poemed_target_unresolvable" %in% seen$classes)
  expect_true("poemed_target_not_reached" %in% seen$classes)
})

test_that("a target no candidate reaches gives NA, a classed warning, and a printed verdict", {
  # The grid `big` selects nothing, so no replication rejects.
  expect_warning(
    plan <- ss_power_pe_mediation(p = 10, n_grid = c(40, 60), n_rep = 5,
                                  lambda_grid = big, seed = 113),
    class = "poemed_target_not_reached")
  expect_true(is.na(attr(plan, "recommended_n")))
  expect_s3_class(plan, "poemed_ss_plan")
  expect_s3_class(plan, "poemed_tbl")
  w <- tryCatch(ss_power_pe_mediation(p = 10, n_grid = c(40, 60), n_rep = 5,
                                      lambda_grid = big, seed = 113),
                poemed_target_not_reached = function(w) w)
  expect_match(conditionMessage(w),
               "Target power 0.8 .*not reached in the grid; the largest estimated power is 0, at n = 40")
  expect_match(conditionMessage(w), "Extend `n_grid` upward")
  out <- capture.output(print(plan))
  expect_true(any(grepl("^Target power 0.8 for the power-enhanced test: not reached in the grid", out)))
  # A plan that reaches the target prints its recommendation.
  reached <- plan
  attr(reached, "recommended_n") <- 60L
  out <- capture.output(print(reached))
  expect_true(any(grepl("recommended n = 60, the smallest n in the grid that reaches it \\(5 replications per n\\)", out)))
  expect_true(any(grepl("^Outcome: continuous", out)))
  # Printing returns the plan unchanged.
  out <- capture.output(printed <- print(reached))
  expect_identical(printed, reached)
})
