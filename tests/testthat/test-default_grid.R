# The default tuning grid, seq(0.05, 10, length.out = 100) for every outcome
# family; the tuning record; and the warning a fit gives when the HBIC choice
# sits at either end of the grid searched.

default_grid <- seq(0.05, 10, length.out = 100)

test_that("every family searches the default grid and records the choice", {
  # The default is written out in every signature, so it can be edited there.
  for (f in list(pe_mediation, pe_mediation_linear, pe_mediation_logistic,
                 pe_mediation_poisson, pe_power_curve, pe_simulation_study,
                 pe_identification_study, ss_power_pe_mediation))
    expect_equal(eval(formals(f)$lambda_grid), default_grid)
  d <- simulate_mediation_data(n = 150, p = 40, outcome = "continuous",
                               pattern = "homogeneous", c1 = 1, seed = 113)
  fit <- quiet_grid(pe_mediation(d$X, d$Y, d$M, outcome = "continuous"))
  tun <- attr(fit, "tuning")
  expect_equal(tun$lambda_grid, default_grid)
  expect_true(tun$lambda_selected %in% tun$lambda_grid)
  # The reduced model of the benchmark test searches the same grid.
  expect_true(tun$lambda_selected_reduced %in% tun$lambda_grid)
  expect_identical(tun$at_lower_end, tun$lambda_selected <= min(tun$lambda_grid))
  expect_identical(tun$at_upper_end, tun$lambda_selected >= max(tun$lambda_grid))
  expect_output(print(fit), "Tuning parameter \\(HBIC\\): lambda = ")
  fit2 <- quiet_grid(pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                                  lambda_grid = c(0.1, 0.2, 0.3)))
  expect_equal(attr(fit2, "tuning")$lambda_grid, c(0.1, 0.2, 0.3))
  expect_true(attr(fit2, "tuning")$lambda_selected_reduced %in% c(0.1, 0.2, 0.3))
  expect_error(pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                            lambda_grid = NULL),
               "`lambda_grid` must be a non-empty numeric vector")

  db <- simulate_mediation_data(n = 200, p = 30, outcome = "binary",
                                pattern = "homogeneous", c1 = 1, c2 = 1,
                                seed = 113)
  fb <- quiet_grid(pe_mediation(db$X, db$Y, db$M, outcome = "binary"))
  expect_equal(attr(fb, "tuning")$lambda_grid, default_grid)
  expect_null(attr(fb, "tuning")$lambda_selected_reduced)
  dc <- simulate_mediation_data(n = 200, p = 30, outcome = "count",
                                pattern = "homogeneous", c1 = 0.4, c2 = 0.4,
                                seed = 113)
  fc <- quiet_grid(pe_mediation(dc$X, dc$Y, dc$M, outcome = "count"))
  expect_equal(attr(fc, "tuning")$lambda_grid, default_grid)
  expect_null(attr(fc, "tuning")$lambda_selected_reduced)
})

test_that("a choice at either end of the grid warns, and an interior choice does not", {
  d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                               pattern = "contrasting", c1 = 1, seed = 113)
  # On the default grid HBIC picks the smallest value, 0.05.
  w <- expect_warning(fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous"),
                      class = "poemed_grid_boundary")
  expect_match(conditionMessage(w), "lambda = 0.05, the smallest value of `lambda_grid`",
               fixed = TRUE)
  # Both penalized fits of a continuous outcome sit at that value here.
  expect_match(conditionMessage(w), "the full model and the reduced model", fixed = TRUE)
  expect_match(conditionMessage(w), "extending the grid")
  expect_match(conditionMessage(w), "finer partition")
  expect_true(attr(fit, "tuning")$at_lower_end)
  expect_output(print(fit), "(the grid's lower end)", fixed = TRUE)
  # A finer grid that reaches below 0.05 finds an interior minimum, and the
  # fit is then silent.
  expect_no_warning(finer <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                                          lambda_grid = seq(0.01, 0.2, length.out = 20)))
  tun <- attr(finer, "tuning")
  expect_false(tun$at_lower_end)
  expect_false(tun$at_upper_end)
  expect_gt(tun$lambda_selected, 0.01)
  expect_lt(tun$lambda_selected, 0.2)
  expect_identical(attr(finer, "selected_mediators"), attr(fit, "selected_mediators"))
  # A grid that stops below that minimum puts the choice at its largest value.
  w2 <- expect_warning(top <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                                           lambda_grid = c(0.02, 0.03, 0.04)),
                       class = "poemed_grid_boundary")
  expect_match(conditionMessage(w2), "lambda = 0.04, the largest value of `lambda_grid`",
               fixed = TRUE)
  expect_true(attr(top, "tuning")$at_upper_end)
  expect_output(print(top), "(the grid's upper end)", fixed = TRUE)
})

test_that("the binary and count workers warn at an end of the grid too", {
  db <- simulate_mediation_data(n = 200, p = 60, outcome = "binary",
                                pattern = "contrasting", c1 = 1, c2 = 1,
                                seed = 113)
  expect_warning(pe_mediation(db$X, db$Y, db$M, outcome = "binary"),
                 class = "poemed_grid_boundary")
  dc <- simulate_mediation_data(n = 200, p = 60, outcome = "count",
                                pattern = "contrasting", c1 = 0.5, c2 = 0.4,
                                seed = 113)
  # The count example's choice is interior (about 0.25) on the default grid.
  expect_no_warning(fc <- pe_mediation(dc$X, dc$Y, dc$M, outcome = "count"))
  expect_false(attr(fc, "tuning")$at_lower_end)
  expect_warning(pe_mediation(dc$X, dc$Y, dc$M, outcome = "count",
                              lambda_grid = c(0.3, 0.4, 0.5)),
                 class = "poemed_grid_boundary")
})

test_that("an empty fit warns once, as an empty fit, not as a boundary choice", {
  d <- simulate_mediation_data(n = 100, p = 20, outcome = "continuous",
                               pattern = "contrasting", c1 = 0, seed = 113)
  # A grid far above any useful lambda keeps no mediator at its top.
  cnd <- NULL
  fit <- withCallingHandlers(
    pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                 lambda_grid = c(50, 100)),
    warning = function(w) { cnd <<- c(cnd, class(w)[1L]); invokeRestart("muffleWarning") })
  expect_true(attr(fit, "empty_fit"))
  expect_identical(cnd, "poemed_empty_fit")
})

test_that("the Monte Carlo functions and the WHO analysis do not raise the warning", {
  skip_on_cran()
  expect_no_warning(pe_power_curve(n = 100, p = 20, outcome = "continuous",
                                   pattern = "contrasting", c1_grid = 1,
                                   n_rep = 2, seed = 113))
  expect_no_warning(pe_identification_study(n = 100, p = 20,
                                            outcome = "continuous",
                                            pattern = "contrasting", c1_grid = 1,
                                            n_rep = 2, methods = "Bonferroni",
                                            seed = 113))
  # The global infant-mortality model picks 0.1, the smallest value of the
  # article's grid; the function keeps that grid and says nothing about it.
  expect_no_warning(who <- WHO_mediation_analysis("imr", groupings = "global"))
  des <- WHO_mediation_design("imr")
  hand <- quiet_grid(pe_mediation_linear(scale(des$X), des$Y - mean(des$Y), scale(des$M),
                                         Z = des$Z, scale = FALSE,
                                         lambda_grid = seq(0.1, 2, length.out = 100)))
  expect_true(attr(hand, "tuning")$at_lower_end)
  expect_equal(who$pval_hdmm, hand$value[hand$term == "pval_hdmm"])
})

test_that("a lambda at which every SCAD weight is zero returns the unpenalized fit", {
  # With two or three candidate mediators the smallest default lambda leaves
  # every first-step coefficient above a * lambda, so no column is penalized.
  # glmnet refuses an all-zero penalty factor; the workers return the
  # unpenalized fit instead (least squares through the origin for the linear
  # model), so the fit keeps every candidate and the run does not stop.
  d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                               pattern = "contrasting", c1 = 1, seed = 113)
  fit <- quiet_grid(pe_mediation(d$X, d$Y, d$M[, 1:2], outcome = "continuous"))
  expect_identical(pe_mediators(fit)$mediator, 1:2)
  expect_equal(attr(fit, "tuning")$lambda_selected, 0.05)
  sc <- POEMED:::.pe_scale_data(d$X, d$Y, d$M[, 1:2], NULL, center_y = TRUE)
  lla <- POEMED:::.LLA_h1(sc$X, sc$Y, sc$M, 0.05, 200, 2, 1)
  expect_equal(as.numeric(lla),
               unname(stats::lm.fit(cbind(sc$M, sc$X), sc$Y)$coefficients))
  expect_s3_class(quiet_grid(pe_mediation(d$X, d$Y, d$M[, 1:3],
                                          outcome = "continuous")), "poemed_tbl")
  dc <- simulate_mediation_data(n = 200, p = 60, outcome = "count",
                                pattern = "contrasting", c1 = 0.5, c2 = 0.4,
                                seed = 113)
  fc <- quiet_grid(pe_mediation(dc$X, dc$Y, dc$M[, 1:2], outcome = "count"))
  expect_identical(pe_mediators(fc)$mediator, 1:2)
  expect_s3_class(quiet_grid(pe_mediation(dc$X, dc$Y, dc$M[, 1:2], outcome = "count",
                                          lambda_grid = c(0.001, 0.01))), "poemed_tbl")
  # The full count design is unchanged by the fallback (it never reaches it).
  full <- quiet_grid(pe_mediation(dc$X, dc$Y, dc$M, outcome = "count"))
  expect_identical(as.integer(attr(full, "selected_mediators")), 4:5)
})
