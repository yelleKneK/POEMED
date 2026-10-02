test_that("the power curve returns the documented structure", {
  set.seed(113)
  pc <- pe_power_curve(n = 80, p = 30, outcome = "continuous",
                       pattern = "contrasting", c1_grid = c(0, 0.6),
                       n_rep = 8, lambda_grid = fast_grid)
  expect_s3_class(pc, "poemed_tbl")
  expect_equal(nrow(pc), 2L)
  expect_true(all(c("c1", "rejection_hdmm", "rejection_pe", "n_valid") %in%
                    names(pc)))
  # Rejection rates are proportions in [0, 1].
  expect_true(all(pc$rejection_hdmm >= 0 & pc$rejection_hdmm <= 1))
  expect_true(all(pc$rejection_pe >= 0 & pc$rejection_pe <= 1))
})

test_that("under contrasting mediation the PE test out-rejects the Wald test as c1 grows", {
  skip_on_cran()
  set.seed(113)
  # The article's design (n = 300, p = 500) at c1 = 1. The article's
  # Figure 3(b) draws the PE test rejecting in about 77 percent of
  # replications and the benchmark in about 15. On the article's own grid
  # (20 values from 0.2 to 0.39, with 0.27 to 0.46 for the reduced model,
  # which the package no longer searches separately) the PE rate agrees
  # (0.770 over 1,000
  # replications with seed = 113 and 4 cores), while the benchmark rejects
  # at a rate of 0.275, above the figure, because the penalized fit
  # sometimes keeps only part of the canceling set. This block runs on the
  # default grid, where the gap is wider still: 0.975 against 0.2 in the 40
  # replications below (seed 113), so the margin of 0.3 sits several Monte
  # Carlo standard errors inside it.
  pc <- pe_power_curve(n = 300, p = 500, outcome = "continuous",
                       pattern = "contrasting", c1_grid = c(0, 1),
                       n_rep = 40)
  null_row <- pc[pc$c1 == 0, ]
  alt_row  <- pc[pc$c1 == 1, ]
  # Size: both tests are near nominal under the global null (generous band).
  expect_lt(null_row$rejection_pe, 0.2)
  # Power: the PE test rejects far more often than the Wald test under the
  # contrasting alternative, which is the article's headline result.
  expect_gt(alt_row$rejection_pe, alt_row$rejection_hdmm + 0.3)
})

test_that("the power curve validates alpha_level", {
  expect_error(pe_power_curve(n = 50, p = 20, c1_grid = 0, n_rep = 2,
                              alpha_level = 0), "in \\(0, 1\\)")
})

test_that("the study functions forward the tuning grid and count empty fits", {
  set.seed(113)
  big <- c(5, 10)                                  # selects nothing at any n
  pc <- pe_power_curve(n = 80, p = 30, outcome = "continuous",
                       pattern = "contrasting", c1_grid = c(0, 1), n_rep = 4,
                       lambda_grid = big)
  expect_true("n_empty" %in% names(pc))
  expect_equal(pc$n_empty, c(4L, 4L))
  expect_equal(pc$n_valid, c(4L, 4L))
  expect_equal(pc$rejection_pe, c(0, 0))          # an empty fit never rejects
  expect_equal(attr(pc, "lambda_grid"), big)
  # The default records the package grid.
  pc2 <- pe_power_curve(n = 80, p = 30, outcome = "continuous",
                        pattern = "contrasting", c1_grid = 1, n_rep = 2)
  expect_equal(attr(pc2, "lambda_grid"), seq(0.05, 10, length.out = 100))
  # Two replications cannot resolve the 0.8 target, and an empty grid
  # reaches no target; both classed warnings are expected here.
  ss <- withCallingHandlers(
    ss_power_pe_mediation(outcome = "continuous", pattern = "contrasting",
                          p = 30, n_grid = c(60, 80), n_rep = 2,
                          lambda_grid = big, seed = 113),
    poemed_target_unresolvable = function(w) invokeRestart("muffleWarning"),
    poemed_target_not_reached = function(w) invokeRestart("muffleWarning"))
  expect_equal(ss$n_empty, c(2L, 2L))
  st <- pe_simulation_study(n = 80, p = 30, outcome = "continuous",
                            c1_grid = 1, n_rep = 2, lambda_grid = big,
                            seed = 113)
  expect_true(all(st$n_empty == 2L))
  expect_equal(attr(st, "lambda_grid"), big)
})

test_that("the study functions reject invalid designs instead of returning NaN", {
  expect_error(pe_power_curve(n = 60, p = 20, n_rep = 3, error_level = 2),
               "`error_level` must be")
  expect_error(pe_power_curve(n = 60, p = 20, n_rep = 0), "`n_rep` must be")
  expect_error(pe_power_curve(n = 0, p = 20, n_rep = 3), "`n` must be")
  expect_error(pe_power_curve(n = 60, p = 20, n_rep = 3, c1_grid = c(NA, 1)),
               "`c1_grid` must be")
  expect_error(pe_power_curve(n = 60, p = 20, n_rep = 3,
                              outcome_args = list(seed = 1)),
               "may not contain `seed`")
  expect_error(pe_power_curve(n = 60, p = 20, n_rep = 3,
                              outcome_args = list(n = 40)),
               "may not contain `n`")
  expect_error(pe_power_curve(n = 60, p = 20, n_rep = 3, lambda_grid = c(0, 1)),
               "`lambda_grid` must be")
  expect_error(pe_simulation_study(n = 60, p = 20, n_rep = 3, error_level = 0),
               "`error_level` must be")
  expect_error(ss_power_pe_mediation(p = 20, n_grid = c(50, 60), n_rep = 3,
                                     c1 = c(0.5, 1)),
               "`c1` must be a single")
  expect_error(ss_power_pe_mediation(p = 20, n_grid = c(50, 60), n_rep = 3,
                                     error_level = 2),
               "`error_level` must be")
  expect_error(pe_identification_study(n = 60, p = 20, n_rep = 3,
                                       error_level = 2),
               "`error_level` must be")
})

test_that("every outcome family defaults to the package grid at every design", {
  # The pages say so: 100 values from 0.05 to 10, whatever n, p, and family.
  default <- seq(0.05, 10, length.out = 100)
  expect_equal(eval(formals(pe_power_curve)$lambda_grid), default)
  for (oc in c("continuous", "count")) {
    pc <- pe_power_curve(n = 60, p = 15, outcome = oc, c1_grid = 0,
                         n_rep = 1, seed = 113)
    expect_equal(attr(pc, "lambda_grid"), default, label = oc)
  }
})

test_that("with cores above 1 the grid points share their random draws, as documented", {
  # ?pe_power_curve and ?pe_identification_study say that with cores > 1
  # every point of c1_grid replays the same n_rep draws (common random
  # numbers), while cores = 1 draws fresh data at each point. Two equal grid
  # points make that visible: identical rows in parallel, different rows in
  # serial. If this ever fails, the help pages must change with the code.
  skip_on_cran()
  skip_on_os("windows")
  grid <- seq(0.05, 1, length.out = 20)
  run <- function(cores)
    pe_identification_study(n = 120, p = 30, pattern = "contrasting",
                            c1_grid = c(1, 1), n_rep = 10, methods = "BH",
                            lambda_grid = grid, cores = cores, seed = 113)
  par <- as.data.frame(unclass(run(2L)))
  ser <- as.data.frame(unclass(run(1L)))
  expect_identical(unlist(par[1, -2]), unlist(par[2, -2]))
  expect_false(identical(unlist(ser[1, -2]), unlist(ser[2, -2])))
})
