# The fit must not depend on the direction in which the outcome or a
# mediator happens to be coded. Negating Y negates the total indirect effect
# and leaves every test statistic, p-value, and the selected set unchanged;
# negating a mediator leaves everything unchanged. This failed before the
# SCAD derivative in the local linear approximation was taken at the
# coefficient's magnitude (a coefficient of -2 was penalized as if it were 0).

stats_of <- function(fit) {
  v <- function(t) fit$value[fit$term == t]
  list(stat_hdmm = v("stat_hdmm"), j_pe = v("j_pe"), stat_pe = v("stat_pe"),
       pval_pe = v("pval_pe"), pval_hdmm = v("pval_hdmm"),
       selected = sort(as.integer(attr(fit, "selected_mediators"))))
}

test_that("the SCAD derivative is symmetric in the sign of the coefficient", {
  z <- c(-3, -2, -0.6, -0.1, 0, 0.1, 0.6, 2, 3)
  expect_identical(POEMED:::.deSCAD(z, 0.5), POEMED:::.deSCAD(abs(z), 0.5))
  expect_identical(POEMED:::.deSCAD(-2, 0.5), 0)
})

test_that("continuous outcome: recoding Y or a mediator changes no statistic", {
  set.seed(113)
  d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                               pattern = "contrasting", c1 = 1)
  f0 <- suppressWarnings(pe_mediation(d$X, d$Y, d$M))
  fy <- suppressWarnings(pe_mediation(d$X, -d$Y, d$M))
  M2 <- d$M; flip <- seq(2, ncol(M2), by = 2); M2[, flip] <- -M2[, flip]
  fm <- suppressWarnings(pe_mediation(d$X, d$Y, M2))
  expect_equal(stats_of(fy), stats_of(f0), tolerance = 1e-10)
  expect_equal(stats_of(fm), stats_of(f0), tolerance = 1e-10)
  tie <- function(f) f$value[f$term == "total_indirect_effect"]
  expect_equal(tie(fy), -tie(f0), tolerance = 1e-10)
  expect_equal(tie(fm), tie(f0), tolerance = 1e-10)
  expect_true(length(stats_of(f0)$selected) > 0)
})

test_that("count and binary outcomes: recoding a mediator changes no statistic", {
  skip_on_cran()
  grid <- seq(0.05, 1, length.out = 20)
  for (nm in c("example_count", "example_binary")) {
    d <- get(nm, envir = asNamespace("POEMED"))
    outcome <- if (nm == "example_count") "count" else "binary"
    M2 <- d$M; flip <- seq(1, ncol(M2), by = 2); M2[, flip] <- -M2[, flip]
    f0 <- suppressWarnings(pe_mediation(d$X, d$Y, d$M, Z = d$Z, outcome = outcome, lambda_grid = grid))
    fm <- suppressWarnings(pe_mediation(d$X, d$Y, M2, Z = d$Z, outcome = outcome, lambda_grid = grid))
    expect_equal(stats_of(fm), stats_of(f0), tolerance = 1e-8, label = nm)
  }
})
