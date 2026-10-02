# Numerical anchors: the three shipped example data sets fit at the grid the
# reference implementation searched, pinned to the values the reference
# implementation produced.
#
# Provenance: the review-era PEmediation package (version 0.1.0, the code that
# produced the article's numbers; PEHDGMM_linear / _logistic / _poisson with
# lamb_grid = seq(0.05, 1, length.out = 20)) matched POEMED bitwise on every one
# of these quantities on 2026-09-24 (tools/oracle_checks.R, run live by the
# release gate, plus a 51-case extension covering Z = NULL, BH, BY, and two
# exposures), and the values below were captured from POEMED 1.0.0 on
# 2026-09-25 at full double precision. The tolerance is the oracle script's own
# criterion, 1e-8 relative, for the statistics; 1e-6 relative for the
# p-values (a far upper tail computed with pchisq(lower.tail = FALSE)). A
# change to any estimation internal that moves a number fails here, which the
# shape and identity tests elsewhere cannot catch (nine one-line mutants of the
# internals passed all of them before these pins were added). No interval is
# pinned: the reference implementation reports none and the fit carries none
# (see test-improvements.R).

oracle_grid <- seq(0.05, 1, length.out = 20)
v <- function(tab, t) tab$value[tab$term == t]

pin <- function(fit, stat_hdmm, pval_hdmm, stat_pe, j_pe, pval_pe, beta,
                selected) {
  expect_equal(v(fit, "stat_hdmm"), stat_hdmm, tolerance = 1e-8)
  expect_equal(v(fit, "pval_hdmm"), pval_hdmm, tolerance = 1e-6)
  expect_equal(v(fit, "stat_pe"), stat_pe, tolerance = 1e-8)
  expect_equal(v(fit, "j_pe"), j_pe, tolerance = 1e-8)
  expect_equal(v(fit, "pval_pe"), pval_pe, tolerance = 1e-6)
  expect_equal(v(fit, "total_indirect_effect"), beta, tolerance = 1e-8)
  expect_identical(as.integer(attr(fit, "selected_mediators")), selected)
  expect_equal(v(fit, "n_selected_mediators"), length(selected))
  expect_equal(v(fit, "stat_pe"), v(fit, "stat_hdmm") + v(fit, "j_pe"))
  expect_equal(nrow(fit), 10L)
}

test_that("the continuous example reproduces the reference implementation", {
  d <- example_continuous
  fit <- pe_mediation(d$X, d$Y, d$M, Z = d$Z, outcome = "continuous",
                      lambda_grid = oracle_grid)
  pin(fit,
      stat_hdmm = 4.65332089289092, pval_hdmm = 0.0309935316324783,
      stat_pe = 457.768594595678, j_pe = 453.115273702787,
      pval_pe = 1.47057154809541e-101, beta = 0.302194895933883,
      selected = c(4L, 5L))
  expect_equal(attr(fit, "tuning")$lambda_selected, 0.1)
})

test_that("the binary example reproduces the reference implementation", quiet_grid({
  d <- example_binary
  fit <- pe_mediation(d$X, d$Y, d$M, Z = d$Z, outcome = "binary",
                      lambda_grid = oracle_grid)
  pin(fit,
      stat_hdmm = 8.00887349523881, pval_hdmm = 0.00465486854375112,
      stat_pe = 357.063306636225, j_pe = 349.054433140986,
      pval_pe = 1.22757397946817e-79, beta = 0.402158518540897,
      selected = 5L)
  expect_equal(attr(fit, "tuning")$lambda_selected, 0.05)
}))

test_that("the count example reproduces the reference implementation", {
  d <- example_count
  fit <- pe_mediation(d$X, d$Y, d$M, Z = d$Z, outcome = "count",
                      lambda_grid = oracle_grid)
  pin(fit,
      stat_hdmm = 7.69601289367094, pval_hdmm = 0.00553429429876433,
      stat_pe = 1085.52812390805, j_pe = 1077.83211101438,
      pval_pe = 4.61618579434415e-238, beta = 0.262067430769173,
      selected = 5L)
  expect_equal(attr(fit, "tuning")$lambda_selected, 0.3)
})
