# Regression tests for ss_power_pe_mediation(): the default (unspecified)
# pattern call must work, every supported pattern must be accepted, and the
# RNG state must be restored on exit in both the existing- and absent-seed
# cases (the @param seed contract).

# These small calls use two replications, which cannot resolve the default
# 0.8 target and rarely reach it, so the planner's two classed warnings are
# expected and muffled here; test-study_guards.R makes each one fire.
quiet_plan <- function(expr) {
  withCallingHandlers(expr,
    poemed_target_unresolvable = function(w) invokeRestart("muffleWarning"),
    poemed_target_not_reached = function(w) invokeRestart("muffleWarning"))
}

test_that("ss_power_pe_mediation runs with the pattern left at its default", {
  plan <- quiet_plan(ss_power_pe_mediation(p = 10, n_grid = c(40, 60),
                                           n_rep = 2, seed = 113))
  expect_s3_class(plan, "poemed_tbl")
  expect_s3_class(plan, "poemed_ss_plan")
  expect_equal(nrow(plan), 2L)
  expect_identical(plan$n, c(40L, 60L))
})

test_that("ss_power_pe_mediation accepts each supported pattern", {
  for (pat in c("homogeneous", "contrasting", "heterogeneous"))
    expect_s3_class(
      quiet_plan(ss_power_pe_mediation(pattern = pat, p = 10, n_grid = 40,
                                       n_rep = 2, lambda_grid = fast_grid,
                                       seed = 1)),
      "poemed_tbl")
})

test_that("the candidate sample sizes are reported in increasing order", {
  plan <- quiet_plan(ss_power_pe_mediation(p = 10, n_grid = c(60, 40),
                                           n_rep = 2, seed = 113))
  expect_identical(plan$n, c(40L, 60L))
})

test_that("ss_power_pe_mediation restores the RNG state on exit", {
  # Absent-seed case: .Random.seed must remain absent afterwards.
  if (exists(".Random.seed", envir = globalenv()))
    rm(".Random.seed", envir = globalenv())
  invisible(quiet_plan(ss_power_pe_mediation(p = 10, n_grid = 40, n_rep = 2,
                                             seed = 113)))
  expect_false(exists(".Random.seed", envir = globalenv()))

  # Existing-seed case: the caller's stream is left untouched.
  set.seed(999)
  before <- get(".Random.seed", envir = globalenv())
  invisible(quiet_plan(ss_power_pe_mediation(p = 10, n_grid = 40, n_rep = 2,
                                             seed = 113)))
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

test_that("the help-page example recommends a sample size, as its comments say", {
  # The example on ?ss_power_pe_mediation, run as written: its comments say
  # the printed plan ends with the recommendation and that the recommended
  # n is 150. If this fails, the example's comments must change with it.
  skip_on_cran()
  set.seed(113)
  plan <- ss_power_pe_mediation(outcome = "continuous", pattern = "contrasting",
                                p = 50, n_grid = c(150, 250), n_rep = 5,
                                lambda_grid = seq(0.05, 2, length.out = 20))
  expect_identical(attr(plan, "recommended_n"), 150L)
  out <- capture.output(print(plan))
  expect_match(out[length(out)],
               "^Target power 0.8 for the power-enhanced test: recommended n = 150")
})
