# Tests for the real-data-motivated heterogeneous setting (simulate_guo_mediation
# and the guo_calibration constants).

test_that("guo_calibration coefficients give the documented total indirect effect", {
  expect_equal(length(guo_calibration$alpha_m), 11L)
  expect_equal(length(guo_calibration$Gamma_x), 11L)
  expect_equal(dim(guo_calibration$Gamma_z), c(11L, 9L))
  # beta_per_c1 is exactly sum(Gamma_x * alpha_m) from the shipped (rounded)
  # coefficients, the value the package actually produces.
  expect_equal(guo_calibration$beta_per_c1,
               sum(guo_calibration$Gamma_x * guo_calibration$alpha_m))
  # It prints as -1.597 to four decimals (the headline number a user sees);
  # the article's -1.5977 comes from the unrounded coefficients.
  expect_equal(round(guo_calibration$beta_per_c1, 4), -1.597)
  expect_equal(guo_calibration$beta_per_c1, -1.5977, tolerance = 1e-3)
})

test_that("simulate_guo_mediation beta equals c1 * beta_per_c1 exactly", {
  for (k in c(0.3, 0.5, 1)) {
    d <- simulate_guo_mediation(c1 = k, seed = 1)
    expect_equal(d$beta, k * guo_calibration$beta_per_c1)
  }
})

test_that("simulate_guo_mediation matches the calibrated dimensions and beta", {
  d <- simulate_guo_mediation(c1 = 0.5, seed = 113)
  expect_equal(dim(d$M), c(85L, 1008L))
  expect_equal(dim(d$X), c(85L, 1L))
  expect_null(d$Z)
  expect_equal(d$active_mediators, 1:11)
  expect_equal(d$beta, guo_calibration$beta_per_c1 * 0.5, tolerance = 1e-8)
})

test_that("there is no mediation at the null, and confounders are supported", {
  expect_length(simulate_guo_mediation(c1 = 0, seed = 1)$active_mediators, 0L)
  dc <- simulate_guo_mediation(c1 = 0.5, confounders = TRUE, seed = 7)
  expect_equal(dim(dc$Z), c(85L, 8L))   # 8 covariates returned (intercept absorbed)
  expect_equal(dc$beta, guo_calibration$beta_per_c1 * 0.5, tolerance = 1e-8)
})

test_that("alpha_m variants and a custom vector are accepted", {
  d1 <- simulate_guo_mediation(c1 = 0.5, alpha_m = "contrasting_like", seed = 1)
  expect_equal(d1$alpha_m[1:4], c(1, -0.5, 0.4, -0.3))
  d2 <- simulate_guo_mediation(c1 = 0.5, alpha_m = rep(0.5, 11), seed = 1)
  expect_equal(unique(d2$alpha_m[1:11]), 0.5)
  expect_error(simulate_guo_mediation(c1 = 0.5, alpha_m = 1:3),
               "must have length 11")
})

test_that("the power-enhanced screen fires on the heterogeneous Guo signal", quiet_grid({
  # Fits several models; the release gate and continuous integration run it
  # in full, and CRAN's machines skip it to keep the check short.
  skip_on_cran()
  set.seed(113)
  d <- simulate_guo_mediation(c1 = 1, n = 150)
  f <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                    lambda_grid = seq(0.1, 10, length.out = 50))
  # The screen selects mediators (J_m > 0), so PE strictly improves
  # on the benchmark Wald test under this mixed-sign heterogeneous signal.
  expect_gte(length(attr(f, "selected_mediators")), 1L)
  expect_lt(f$value[f$term == "pval_pe"], f$value[f$term == "pval_hdmm"])
}))

test_that("guo_calibration equals the constants the supplement prints", {
  # Section S.4.5 of the supplement of Yu and Kelley (in press) prints the
  # Guo et al. (2022) estimates to three decimals; Section S.4.4 gives the
  # designed outcome-mediator coefficients and their two alternatives. Each
  # constant is pinned value by value, so a transcription error in any one
  # of them (the confounder constants included) fails here.
  expect_identical(guo_calibration$p, 1008L)
  expect_identical(guo_calibration$s, 9L)
  expect_identical(guo_calibration$n, 85L)
  expect_identical(guo_calibration$locations, 1:11)
  expect_identical(guo_calibration$alpha_m,
                   c(1.0, 0.9, 0.8, -0.9, -0.8, -0.7, 0.6, 0.5, 0.4, 0.3, 0.2))
  expect_identical(guo_calibration$Gamma_x,
                   c(-0.251, -0.221, -0.233, 0.251, 0.295, 0.282,
                     -0.332, -0.359, 0.335, -0.345, 0.234))
  expect_identical(guo_calibration$alpha_z,
                   c(-0.336, -0.070, 0.665, 0.278, 0.315, 0.201, 0.173,
                     0.510, 0.315))
  expect_identical(guo_calibration$alpha_m_estimated,
                   c(0.166, 0.243, 0.248, -0.049, -0.294, -0.187, 0.148,
                     0.087, 0.112, 0.223, 0.165))
  Gamma_z <- matrix(c(
    -0.045,  0.076,  0.089,  0.127, -0.408, -0.233, -0.104, -0.442, -0.242,
    -0.197,  0.100,  0.390, -0.357, -0.310, -0.195, -0.270, -0.497, -0.419,
    -0.076,  0.149,  0.151, -0.590, -0.813, -0.242, -0.273, -1.217, -0.614,
    -0.052,  0.048,  0.103,  0.115,  0.065, -0.065,  0.114,  0.199,  0.003,
    -0.033, -0.184,  0.065,  0.156,  0.008,  0.116,  0.070,  0.365,  0.287,
    -0.012, -0.095,  0.023,  0.112,  0.200, -0.147,  0.105,  0.240, -0.080,
     0.203,  0.396, -0.400,  0.377,  0.674,  0.543,  0.452,  1.180,  0.671,
     0.017,  0.234, -0.034,  0.167, -0.411, -0.063,  0.039, -0.260, -0.091,
    -0.364,  0.082,  0.719, -0.252, -0.517,  0.007, -0.066, -0.793, -0.349,
    -0.001,  0.098,  0.001,  0.040, -0.052,  0.064, -0.040, -0.076, -0.070,
    -0.099,  0.022,  0.196, -0.103,  0.468,  0.138,  0.213,  0.491,  0.122),
    nrow = 11, ncol = 9, byrow = TRUE)
  expect_identical(unname(guo_calibration$Gamma_z), Gamma_z)
  expect_identical(guo_calibration$alpha_m_variants$homogeneous_like,
                   c(1, 0.8, 0.6, 0.4, 0.2, rep(0, 6)))
  expect_identical(guo_calibration$alpha_m_variants$contrasting_like,
                   c(1, -0.5, 0.4, -0.3, rep(0, 7)))
})

test_that("active_mediators lists the loci nonzero in both paths", {
  # All eleven under the default; the first five and first four under the
  # two alternative sets, whose trailing loci have no outcome effect; only
  # the loci with a nonzero entry in a numeric alpha_m; none at c1 = 0.
  expect_identical(simulate_guo_mediation(c1 = 1, n = 20, seed = 113)$active_mediators,
                   1:11)
  expect_identical(simulate_guo_mediation(c1 = 1, n = 20, seed = 113,
                                          alpha_m = "homogeneous_like")$active_mediators,
                   1:5)
  expect_identical(simulate_guo_mediation(c1 = -0.5, n = 20, seed = 113,
                                          alpha_m = "contrasting_like")$active_mediators,
                   1:4)
  d <- simulate_guo_mediation(c1 = 1, n = 20, seed = 113,
                              alpha_m = c(1, 0, 0.5, rep(0, 8)))
  expect_identical(d$active_mediators, c(1L, 3L))
  expect_identical(d$active_mediators,
                   which(d$alpha_m != 0 & as.numeric(d$Gamma_x) != 0))
  expect_length(simulate_guo_mediation(c1 = 0, n = 20, seed = 113,
                                       alpha_m = "homogeneous_like")$active_mediators,
                0L)
})

test_that("simulate_guo_mediation validates its arguments by name", {
  expect_error(simulate_guo_mediation(c1 = NA), "`c1` must be a single finite number")
  expect_error(simulate_guo_mediation(c1 = Inf), "`c1` must be")
  expect_error(simulate_guo_mediation(c1 = c(0.5, 1)), "`c1` must be")
  expect_error(simulate_guo_mediation(c1 = "a"), "`c1` must be")
  expect_error(simulate_guo_mediation(1, c2 = NA), "`c2` must be a single finite number")
  expect_error(simulate_guo_mediation(1, n = c(10, 20)),
               "`n` must be a single whole number of at least 2")
  expect_error(simulate_guo_mediation(1, n = 0), "`n` must be")
  expect_error(simulate_guo_mediation(1, n = 2.5), "`n` must be")
  for (v in list(NA, "yes", c(TRUE, FALSE)))
    expect_error(simulate_guo_mediation(1, confounders = v),
                 "`confounders` must be TRUE or FALSE")
  expect_error(simulate_guo_mediation(1, alpha_m = rep(NA_real_, 11)),
               "A numeric `alpha_m` must hold finite values")
  expect_error(simulate_guo_mediation(1, alpha_m = "setting1"),
               "`alpha_m` must be one of \"setting0\", \"homogeneous_like\", \"contrasting_like\"")
  expect_error(simulate_guo_mediation(1, alpha_m = c("setting0", "homogeneous_like")),
               "`alpha_m` must be one of")
})
