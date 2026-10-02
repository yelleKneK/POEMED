test_that("the simulator returns the documented structure and truth", {
  set.seed(113)
  d <- simulate_mediation_data(n = 100, p = 40, outcome = "continuous",
                               pattern = "homogeneous", d = 2)
  expect_equal(dim(d$X), c(100L, 1L))
  expect_equal(dim(d$M), c(100L, 40L))
  expect_length(d$Y, 100L)
  expect_equal(dim(d$Z), c(100L, 2L))
  expect_length(d$alpha_m, 40L)
  # beta is the total indirect effect Gamma_x %*% alpha_m, and the active set
  # is exactly the mediators nonzero on both paths.
  expect_equal(d$beta, as.numeric(d$Gamma_x %*% d$alpha_m))
  expect_equal(d$active_mediators,
               which(d$alpha_m != 0 & apply(d$Gamma_x != 0, 2, any)))
})

test_that("the contrasting pattern has mixed-sign active mediators", {
  for (oc in c("continuous", "binary", "count")) {
    set.seed(113)
    d <- simulate_mediation_data(n = 80, p = 30, outcome = oc,
                                 pattern = "contrasting")
    expect_gt(length(d$active_mediators), 1L)
    # Individual indirect effects (path product) of the active mediators carry
    # both signs: this is what "contrasting" means. The continuous and binary
    # presets cancel to an exactly zero total with the default tau; the count
    # preset (from the article) is mixed-sign but does not exactly cancel.
    ind <- (d$Gamma_x %*% diag(d$alpha_m))[1, d$active_mediators]
    expect_true(any(ind > 0) && any(ind < 0), label = paste("signs for", oc))
  }
})

test_that("the continuous and binary contrasting presets cancel to zero total", {
  for (oc in c("continuous", "binary")) {
    set.seed(113)
    d <- simulate_mediation_data(n = 80, p = 30, outcome = oc,
                                 pattern = "contrasting")
    expect_lt(abs(d$beta), 1e-8, label = paste("beta for", oc))
  }
})

test_that("outcome types produce correctly typed responses", {
  set.seed(113)
  b <- simulate_mediation_data(n = 80, p = 30, outcome = "binary",
                               pattern = "homogeneous")
  expect_true(all(b$Y %in% c(0, 1)))
  set.seed(113)
  k <- simulate_mediation_data(n = 80, p = 30, outcome = "count",
                               pattern = "homogeneous")
  expect_true(all(k$Y >= 0 & k$Y == round(k$Y)))
})

test_that("c1 = 0 yields the global null (no exposure-mediator signal)", {
  set.seed(113)
  d <- simulate_mediation_data(n = 80, p = 30, outcome = "continuous",
                               pattern = "homogeneous", c1 = 0)
  expect_true(all(d$Gamma_x == 0))
  expect_equal(d$beta, 0)
})

test_that("seeding is reproducible and restores the caller's RNG state", {
  d1 <- simulate_mediation_data(n = 50, p = 20, outcome = "continuous",
                                pattern = "homogeneous", seed = 7)
  d2 <- simulate_mediation_data(n = 50, p = 20, outcome = "continuous",
                                pattern = "homogeneous", seed = 7)
  expect_identical(d1$M, d2$M)
  expect_identical(d1$Y, d2$Y)
  # The caller's stream is untouched by a seeded call.
  set.seed(42); before <- .Random.seed
  invisible(simulate_mediation_data(n = 30, p = 10, outcome = "continuous",
                                    pattern = "homogeneous", seed = 99))
  expect_identical(.Random.seed, before)
})

test_that("a user-supplied alpha_m overrides the preset", {
  set.seed(113)
  am <- c(2, -2, rep(0, 28))
  d <- simulate_mediation_data(n = 60, p = 30, outcome = "continuous",
                               pattern = "homogeneous", alpha_m = am)
  expect_identical(d$alpha_m, am)
  expect_error(simulate_mediation_data(n = 60, p = 30, alpha_m = c(1, 2)),
               "length")
})

test_that("the simulator validates its design arguments by name", {
  expect_error(simulate_mediation_data(n = 0, p = 20), "`n` must be")
  expect_error(simulate_mediation_data(n = 50, p = 1.5), "`p` must be")
  expect_error(simulate_mediation_data(n = 50, p = 20, c1 = "a"), "`c1` must be")
  expect_error(simulate_mediation_data(n = 50, p = 20, c2 = c(1, 2)), "`c2` must be")
  expect_error(simulate_mediation_data(n = 50, p = 20, sigma_y = 0), "`sigma_y` must be")
  expect_error(simulate_mediation_data(n = 50, p = 20, sigma_y = Inf), "`sigma_y` must be")
  expect_error(simulate_mediation_data(n = 50, p = 20, d = -1), "`d` must be")
  expect_error(simulate_mediation_data(n = 50, p = 20, q = 0), "`q` must be")
  expect_error(simulate_mediation_data(n = 50, p = 20, rho = 1), "`rho` must")
  expect_error(simulate_mediation_data(n = 50, p = 20, rho = NA), "`rho` must")
})

test_that("the six alpha_m presets are the article's", {
  # Section 4.1 of Yu and Kelley (in press) for the continuous outcome and
  # Section 4.2 for the binary (logistic) and count (Poisson) outcomes:
  # setting (i) is the homogeneous preset and setting (ii) the contrasting
  # one. Each preset is padded with zeros to length p.
  presets <- list(
    "continuous/homogeneous" = c(1, 0.8, 0.6, 0.4, 0.2),
    "continuous/contrasting" = c(1, -0.5, 0.4, -0.3),
    "binary/homogeneous"     = c(3, 1.5, 0, 0, 2),
    "binary/contrasting"     = c(3, -1.5),
    "count/homogeneous"      = c(0.9, 0.8, 0, 0, 0.7),
    "count/contrasting"      = c(0, 0, 0, 0.8, -0.7))
  for (k in names(presets)) {
    s <- strsplit(k, "/")[[1]]
    a <- simulate_mediation_data(n = 20, p = 10, outcome = s[1],
                                 pattern = s[2], seed = 113)$alpha_m
    nz <- presets[[k]]
    expect_identical(a[seq_along(nz)], nz, label = k)
    expect_true(all(a[-seq_along(nz)] == 0), label = k)
    expect_length(a, 10L)
  }
})

test_that("the contrasting total indirect effect is zero, or -0.03 c1 for the count preset", {
  # ?simulate_mediation_data (the `pattern` argument) says: with the default
  # tau the continuous and binary contrasting presets cancel, while the count
  # preset (0, 0, 0, 0.8, -0.7) gives 0.4 * 0.8 - 0.5 * 0.7 = -0.03 per c1.
  for (oc in c("continuous", "binary")) {
    d <- simulate_mediation_data(n = 40, p = 20, outcome = oc,
                                 pattern = "contrasting", seed = 113)
    expect_equal(d$beta, 0, label = paste("beta for", oc))
  }
  for (c1 in c(1, 0.5)) {
    d <- simulate_mediation_data(n = 40, p = 20, outcome = "count",
                                 pattern = "contrasting", c1 = c1, seed = 113)
    expect_equal(d$beta, -0.03 * c1)
  }
})

test_that("alpha_m, tau, alpha_z, and Gamma_z are validated by name", {
  expect_error(simulate_mediation_data(n = 50, p = 12,
                                       alpha_m = c(NA, rep(0, 11))),
               "`alpha_m` must be a numeric vector of finite values")
  expect_error(simulate_mediation_data(n = 50, p = 12,
                                       alpha_m = c(Inf, rep(0, 11))),
               "`alpha_m` must be")
  expect_error(simulate_mediation_data(n = 50, p = 12,
                                       alpha_m = rep("a", 12)),
               "`alpha_m` must be")
  expect_error(simulate_mediation_data(n = 50, p = 12, alpha_m = rep(0, 11)),
               "`alpha_m` must have length `p` \\(12\\)")
  expect_error(simulate_mediation_data(n = 50, p = 10, tau = c(NA, 1:9)),
               "`tau` must be a 1-by-10 matrix of finite numbers")
  expect_error(simulate_mediation_data(n = 50, p = 10, tau = rep("a", 10)),
               "`tau` must be")
  expect_error(simulate_mediation_data(n = 50, p = 10, tau = 1:9),
               "`tau` must be a 1-by-10 matrix")
  # With two exposures a p-by-q matrix would be silently rearranged.
  expect_error(simulate_mediation_data(n = 50, p = 10, q = 2,
                                       tau = matrix(0.1, 10, 2)),
               "`tau` must be a 2-by-10 matrix")
  expect_error(simulate_mediation_data(n = 50, p = 10, q = 2),
               "For `q > 1` you must supply `tau`")
  expect_error(simulate_mediation_data(n = 50, p = 4, pattern = "contrasting"),
               "`p` must be at least 5 for the default `tau`")
  expect_error(simulate_mediation_data(n = 50, p = 12, d = 2, alpha_z = 1:3),
               "`alpha_z` must be a numeric vector of 2 finite values")
  expect_error(simulate_mediation_data(n = 50, p = 12, d = 2,
                                       alpha_z = c(1, NA)),
               "`alpha_z` must be")
  expect_error(simulate_mediation_data(n = 50, p = 10, d = 2,
                                       Gamma_z = matrix(1, 3, 10)),
               "`Gamma_z` must be a 2-by-10 matrix")
  expect_error(simulate_mediation_data(n = 50, p = 10, d = 2,
                                       Gamma_z = rep(1, 20)),
               "`Gamma_z` must be a 2-by-10 matrix")
  expect_error(simulate_mediation_data(n = 50, p = 10, d = 1,
                                       Gamma_z = c(1, rep(NaN, 9))),
               "`Gamma_z` must be a 1-by-10 matrix of finite numbers")
})

test_that("confounder coefficients supplied with d = 0 stop instead of being ignored", {
  expect_error(simulate_mediation_data(n = 50, p = 12, alpha_z = c(1, 1)),
               "`alpha_z` describes confounders, but `d = 0` asks for none")
  expect_error(simulate_mediation_data(n = 50, p = 12,
                                       Gamma_z = matrix(1, 2, 12)),
               "`Gamma_z` describes confounders")
  expect_error(simulate_mediation_data(n = 50, p = 12, alpha_z = c(5, 5),
                                       Gamma_z = matrix(1, 2, 12)),
               "`alpha_z` and `Gamma_z` describe confounders.*Set `d`")
  # Empty coefficients are consistent with d = 0 and change nothing.
  a <- simulate_mediation_data(n = 30, p = 12, seed = 113)
  b <- simulate_mediation_data(n = 30, p = 12, alpha_z = numeric(0),
                               seed = 113)
  expect_identical(a, b)
})

test_that("valid coefficient shapes are accepted as before", {
  # One exposure: tau as a vector, a row, or a column holding p values.
  tau <- seq(0.05, 1, length.out = 10)
  a <- simulate_mediation_data(n = 30, p = 10, tau = tau, seed = 113)
  b <- simulate_mediation_data(n = 30, p = 10, tau = matrix(tau, 1, 10),
                               seed = 113)
  cc <- simulate_mediation_data(n = 30, p = 10, tau = matrix(tau, 10, 1),
                                seed = 113)
  expect_identical(a, b)
  expect_identical(a, cc)
  # Two exposures: a 2-by-p matrix, or its values in column order.
  tau2 <- matrix(seq(0.05, 1, length.out = 20), 2, 10)
  a <- simulate_mediation_data(n = 30, p = 10, q = 2, tau = tau2, seed = 113)
  b <- simulate_mediation_data(n = 30, p = 10, q = 2, tau = as.vector(tau2),
                               seed = 113)
  expect_identical(a, b)
  expect_identical(a$Gamma_x, tau2)
  # One confounder: Gamma_z as a length-p vector or a 1-by-p matrix.
  g <- seq(0, 1, length.out = 10)
  a <- simulate_mediation_data(n = 30, p = 10, d = 1, alpha_z = 0.4,
                               Gamma_z = g, seed = 113)
  b <- simulate_mediation_data(n = 30, p = 10, d = 1, alpha_z = 0.4,
                               Gamma_z = matrix(g, 1, 10), seed = 113)
  expect_identical(a, b)
  expect_equal(dim(a$Z), c(30L, 1L))
})

test_that("outcome and pattern name themselves when they do not match", {
  expect_error(simulate_mediation_data(n = 50, p = 12, outcome = "poisson"),
               "`outcome` must be one of \"continuous\", \"binary\", \"count\"")
  expect_error(simulate_mediation_data(n = 50, p = 12, pattern = "mixed"),
               "`pattern` must be one of")
  expect_error(simulate_mediation_data(n = 50, p = 12,
                                       outcome = c("binary", "count")),
               "`outcome` must be one of")
})

test_that("the design checker keeps the simulator's arguments and defaults", {
  # The study functions validate a design with .pe_sim_design() before any
  # replication, so its arguments and defaults must stay those of
  # simulate_mediation_data() (without `seed`).
  sim <- formals(simulate_mediation_data)
  chk <- formals(POEMED:::.pe_sim_design)
  expect_identical(names(chk), setdiff(names(sim), "seed"))
  expect_identical(as.list(chk), as.list(sim)[names(chk)])
})
