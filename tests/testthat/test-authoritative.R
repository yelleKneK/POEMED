# Authoritative values: each result checked against a number computed outside
# the code that produced it.
#
# The shape, range, and identity tests elsewhere show that a result has the
# right columns and plausible values; they cannot tell a correct number from a
# wrong one of the same shape. Each test here sets a returned number against an
# independent authority:
#
#   * the Monte Carlo studies (pe_power_curve(), pe_simulation_study(),
#     ss_power_pe_mediation(), pe_identification_study()) against a loop over
#     pe_mediation() written out below, drawn from the same seed and scored by
#     the definitions on their help pages;
#   * pe_mediators(), tidy(), and glance() against the fit's own rows and
#     against least squares fits written out with lm();
#   * the BH and BY screens against the Benjamini-Hochberg and
#     Benjamini-Yekutieli step-up rules written out by hand;
#   * the p-values with two exposures against the chi-square upper tail on
#     q = 2 degrees of freedom, which is exp(-x / 2) in closed form;
#   * pe_mediate() with confounders against pe_mediation() given the
#     confounder matrix built with model.matrix().
#
# Provenance of the pinned numbers. They were computed with POEMED 1.0.0 on
# 2026-09-29 at full double precision. Every pinned statistic, selected set, and
# total indirect effect was also produced by the review-era reference
# implementation (PEmediation 0.1.0: PEHDGMM_linear, PEHDGMM_logistic, and
# PEHDGMM_poisson) given the same data and the same tuning grids, bit for bit
# or to at least 14 significant digits. For a continuous fit on the then
# default grid, the two grids then recorded in the fit's "tuning" attribute
# were passed to the reference as lamb_grid and lamb_grid0; since 2026-10-01
# the record holds one grid, passed as both (see pinned_grid() below). The reference returns only the PE statistic, so
# its Wald statistic was read from a run whose screening level was so small
# that no mediator passed, which leaves the PE statistic equal to the Wald
# statistic. The reference reports its p-value as 1 - pchisq(), which is 0 in
# the far tail, so far-tail p-values are checked against pchisq() with
# lower.tail = FALSE on the log scale instead. Tolerances follow
# test-oracle_pins.R: 1e-8 relative for statistics and effects. The
# continuous fits pinned on what was then the default grid pass that grid
# explicitly (pinned_grid() below), so the pins and their cross-check stand
# under the current default, seq(0.05, 10, length.out = 100).
#
# The Monte Carlo tests need no pins, since their authority is the loop
# written out here. For the record, the PE rejection rates and the
# identification metrics of those designs also agree exactly with the
# reference run on the same simulated data sets.

v <- function(tab, t) tab$value[tab$term == t]
oracle_grid <- seq(0.05, 1, length.out = 20)

# The grid under which the continuous-outcome numbers below were pinned and
# cross-checked against the reference implementation: the article's linear
# grid (0.2 to 0.39) scaled by sqrt(log(p) / n) relative to the article's
# n = 300, p = 500, the package's default for a continuous outcome when the
# pins were first computed. The reduced model of the benchmark test searches
# the same grid (the reference's default lamb_grid0 = lamb_grid). On
# 2026-10-01 every continuous pin was recomputed with that single grid and
# checked again against the reference, which gave the same values to all
# printed digits: the reduced model selects the same set on either grid.
pinned_grid <- function(n, p) {
  ends <- c(0.2, 0.39)
  ref_rate <- sqrt(log(500) / 300)
  rate <- sqrt(log(p) / n)
  seq(ends[1L] / ref_rate * rate, ends[2L] / ref_rate * rate, length.out = 20)
}

# Fit without the classed empty-fit warning (the Monte Carlo studies count
# such fits rather than warn) and without the grid-boundary warning (the
# studies search a fixed grid by design); every other warning still surfaces.
fit_quietly <- function(...) {
  withCallingHandlers(pe_mediation(...),
                      poemed_empty_fit = function(w)
                        invokeRestart("muffleWarning"),
                      poemed_grid_boundary = function(w)
                        invokeRestart("muffleWarning"))
}

# One replication of the power study written out by hand: fit pe_mediation()
# to a simulated data set and record whether each test rejects at alpha_level
# and whether the penalized selection was empty.
hand_rejections <- function(dat, outcome, alpha_level = 0.05, ...) {
  fit <- fit_quietly(dat$X, dat$Y, dat$M, Z = dat$Z, outcome = outcome, ...)
  c(hdmm = v(fit, "pval_hdmm") <= alpha_level,
    pe = v(fit, "pval_pe") <= alpha_level,
    empty = isTRUE(attr(fit, "empty_fit")))
}

# A power curve written out by hand. For each c1, n_rep replications are drawn
# in turn from the current random stream, each one simulated and then fit.
# The rates are the proportions of replications that reject.
hand_curve <- function(n, p, outcome, pattern, c1_grid, n_rep, c2 = 0.5, ...) {
  rows <- lapply(c1_grid, function(c1) {
    r <- vapply(seq_len(n_rep), function(i) {
      dat <- simulate_mediation_data(n = n, p = p, outcome = outcome,
                                     pattern = pattern, c1 = c1, c2 = c2)
      hand_rejections(dat, outcome, ...)
    }, logical(3))
    data.frame(c1 = c1, rejection_hdmm = mean(r["hdmm", ]),
               rejection_pe = mean(r["pe", ]), n_valid = as.integer(n_rep),
               n_empty = sum(r["empty", ]))
  })
  do.call(rbind, rows)
}

expect_same_curve <- function(got, hand) {
  expect_equal(got$c1, hand$c1)
  expect_equal(got$rejection_hdmm, hand$rejection_hdmm)
  expect_equal(got$rejection_pe, hand$rejection_pe)
  expect_identical(got$n_valid, hand$n_valid)
  expect_identical(got$n_empty, hand$n_empty)
}

# The Benjamini-Hochberg step-up rule, and with by = TRUE the
# Benjamini-Yekutieli rule (the level divided by the harmonic sum
# 1 + 1/2 + ... + 1/m): reject the hypotheses with the i smallest p-values,
# where i is the largest rank whose sorted p-value is at most i / m times the
# level. Returns a logical vector in the order of `p`.
step_up <- function(p, level, by = FALSE) {
  m <- length(p)
  if (by) level <- level / sum(1 / seq_len(m))
  o <- order(p)
  passing <- which(p[o] <= seq_len(m) / m * level)
  keep <- logical(m)
  if (length(passing)) keep[o[seq_len(max(passing))]] <- TRUE
  keep
}

# The two-exposure design shared by the q = 2 tests: the loadings of the
# second exposure reverse those of the first.
two_exposure_data <- function(outcome, c1) {
  tau <- rbind(c(0.1, 0.2, 0.3, 0.4, 0.5, rep(0, 25)),
               c(0.5, 0.4, 0.3, 0.2, 0.1, rep(0, 25)))
  simulate_mediation_data(n = 200, p = 30, outcome = outcome,
                          pattern = "contrasting", c1 = c1, q = 2, tau = tau,
                          seed = 113)
}

fit_example <- function(method = "Bonferroni", ...) {
  d <- example_continuous
  pe_mediation(d$X, d$Y, d$M, Z = d$Z, outcome = "continuous",
               method = method, lambda_grid = oracle_grid, ...)
}

# ---------------------------------------------------------------------------
# The Monte Carlo studies against a hand-written loop over pe_mediation().
# ---------------------------------------------------------------------------

test_that("pe_power_curve() rates equal a hand-written loop over pe_mediation()", quiet_grid({
  pc <- pe_power_curve(n = 200, p = 60, outcome = "continuous",
                       pattern = "contrasting", c1_grid = c(1, 1.5),
                       n_rep = 3, lambda_grid = fast_grid, seed = 113)
  set.seed(113)
  hand <- hand_curve(200, 60, "continuous", "contrasting", c(1, 1.5), 3,
                     lambda_grid = fast_grid)
  expect_same_curve(pc, hand)
  # The design separates the two tests, so exchanging the benchmark and
  # power-enhanced columns would fail above.
  expect_false(identical(hand$rejection_hdmm, hand$rejection_pe))
}))

test_that("pe_simulation_study() rows equal hand-written loops, pattern by pattern", quiet_grid({
  skip_on_cran()
  st <- pe_simulation_study(n = 150, p = 30, outcome = "continuous",
                            c1_grid = 1, n_rep = 6, seed = 113)
  # The study seeds once and runs the patterns in turn on one stream.
  set.seed(113)
  hom <- hand_curve(150, 30, "continuous", "homogeneous", 1, 6)
  con <- hand_curve(150, 30, "continuous", "contrasting", 1, 6)
  expect_identical(st$pattern, c("homogeneous", "contrasting"))
  expect_same_curve(st[st$pattern == "homogeneous", ], hom)
  expect_same_curve(st[st$pattern == "contrasting", ], con)
  # The patterns give different rates here, so a row carrying the wrong
  # pattern label would fail above.
  expect_gt(hom$rejection_hdmm, con$rejection_hdmm)
}))

test_that("ss_power_pe_mediation() powers equal a hand-written loop and recommend the smallest n", quiet_grid({
  skip_on_cran()
  plan <- ss_power_pe_mediation(outcome = "continuous",
                                pattern = "contrasting", p = 30,
                                n_grid = c(150, 250, 350), c1 = 1,
                                target_power = 0.5, n_rep = 8, seed = 113)
  set.seed(113)
  hand <- do.call(rbind, lapply(c(150, 250, 350), function(nn)
    hand_curve(nn, 30, "continuous", "contrasting", 1, 8)))
  expect_equal(plan$n, c(150L, 250L, 350L))
  expect_equal(plan$power_pe, hand$rejection_pe)
  expect_equal(plan$power_hdmm, hand$rejection_hdmm)
  expect_identical(plan$n_valid, hand$n_valid)
  expect_identical(plan$n_empty, hand$n_empty)
  # The recommendation is the smallest candidate n whose PE power reaches
  # the target. More than one candidate reaches it here, so the smallest and
  # the largest differ.
  reached <- plan$n[hand$rejection_pe >= 0.5]
  expect_gt(length(reached), 1L)
  expect_identical(attr(plan, "recommended_n"), min(reached))
  expect_identical(attr(plan, "recommended_n"), 150L)
}))

test_that("pe_identification_study() metrics equal their definitions scored by hand", quiet_grid({
  skip_on_cran()
  # Mediator 1 keeps its outcome effect but loses its exposure path (its tau
  # is 0), so it is a true mediator under truth = "outcome_effect" and a false
  # one under truth = "mediation". Its exposure-path p-value is then noise,
  # and the loose error_level lets the screen admit it now and then, so some
  # replications select true and false mediators together.
  tau <- matrix(c(0, 0.2, 0.3, 0.4, 0.5, rep(0, 25)), nrow = 1)
  design <- list(n = 200, p = 30, outcome = "continuous",
                 pattern = "contrasting", c1_grid = c(0, 1), n_rep = 6,
                 error_level = 0.9, outcome_args = list(tau = tau),
                 lambda_grid = pinned_grid(200, 30))
  methods <- c("Bonferroni", "BH", "BY")
  id_med <- do.call(pe_identification_study,
                    c(design, list(truth = "mediation", seed = 113)))
  id_out <- do.call(pe_identification_study,
                    c(design, list(truth = "outcome_effect", seed = 113)))

  # The same replications by hand: simulate, fit once per method, and keep
  # each method's selected set with both definitions of the truth set.
  set.seed(113)
  reps <- lapply(design$c1_grid, function(c1) lapply(seq_len(design$n_rep), function(r) {
    dat <- simulate_mediation_data(n = design$n, p = design$p,
                                   outcome = design$outcome,
                                   pattern = design$pattern, c1 = c1,
                                   tau = tau)
    fits <- lapply(methods, function(m)
      fit_quietly(dat$X, dat$Y, dat$M, outcome = design$outcome, method = m,
                  error_level = design$error_level,
                  lambda_grid = design$lambda_grid))
    list(sets = stats::setNames(lapply(fits, attr, "selected_mediators"), methods),
         empty = isTRUE(attr(fits[[1]], "empty_fit")),
         truth = list(mediation = dat$active_mediators,
                      outcome_effect = which(dat$alpha_m != 0)))
  }))

  # The help page's definitions, per replication: a familywise error is any
  # selected mediator outside the truth set; the false discovery proportion
  # and precision are shares of the selected set (0 when nothing is
  # selected); recall is the share of the truth set selected (undefined when
  # the truth set is empty). The study reports their means.
  score <- function(sel, truth) {
    fp <- length(setdiff(sel, truth)); tp <- length(intersect(sel, truth))
    c(fwer = as.numeric(fp > 0),
      fdr = if (length(sel)) fp / length(sel) else 0,
      precision = if (length(sel)) tp / length(sel) else 0,
      recall = if (length(truth)) tp / length(truth) else NA_real_)
  }
  mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  by_hand <- function(truth_name) {
    do.call(rbind, lapply(seq_along(design$c1_grid), function(g)
      do.call(rbind, lapply(methods, function(m) {
        s <- vapply(reps[[g]], function(x) score(x$sets[[m]], x$truth[[truth_name]]),
                    numeric(4))
        data.frame(c1 = design$c1_grid[g], method = m,
                   fwer = mean_or_na(s["fwer", ]), fdr = mean_or_na(s["fdr", ]),
                   precision = mean_or_na(s["precision", ]),
                   recall = mean_or_na(s["recall", ]),
                   n_valid = length(reps[[g]]),
                   n_empty = sum(vapply(reps[[g]], `[[`, logical(1), "empty")))
      }))))
  }

  for (tr in c("mediation", "outcome_effect")) {
    got <- if (tr == "mediation") id_med else id_out
    hand <- by_hand(tr)
    expect_equal(got$c1, hand$c1)
    expect_identical(got$method, hand$method)
    for (col in c("fwer", "fdr", "precision", "recall"))
      expect_equal(got[[col]], hand[[col]], info = paste(tr, col))
    expect_identical(got$n_valid, hand$n_valid)
    expect_identical(got$n_empty, hand$n_empty)
  }

  # The design exercises every definition: familywise errors under the null
  # with recall undefined there, a false discovery rate below the familywise
  # error rate (a replication mixing true and false selections), precision
  # strictly between 0 and 1, recall that differs between methods, and the
  # two truth definitions disagreeing at c1 = 1.
  med <- by_hand("mediation")
  out <- by_hand("outcome_effect")
  expect_true(all(med$fwer[med$c1 == 0] > 0))
  expect_true(all(is.na(med$recall[med$c1 == 0])))
  expect_true(any(med$fdr > 0 & med$fdr < med$fwer))
  expect_true(any(med$precision > 0 & med$precision < 1))
  expect_gt(length(unique(med$recall[med$c1 == 1])), 1L)
  expect_false(identical(med$precision[med$c1 == 1], out$precision[out$c1 == 1]))
}))

# ---------------------------------------------------------------------------
# pe_mediators(), tidy(), and glance() against the fit and against lm().
# ---------------------------------------------------------------------------

test_that("pe_mediators() and tidy() equal least squares fits written out with lm()", {
  d <- example_continuous
  fit <- fit_example()
  tab <- pe_mediators(fit)
  A <- tab$mediator
  expect_identical(A, 1:5)
  # A t statistic does not change when the package centers and scales the
  # variables, so the raw data serve. The outcome path is the regression of Y
  # on the selected mediators, the exposure, and the confounders; the
  # exposure path is the regression of each mediator on the exposure and the
  # confounders.
  ols <- summary(lm(d$Y ~ d$M[, A] + d$X + d$Z))$coefficients
  t_out <- unname(ols[1 + seq_along(A), "t value"])
  t_exp <- vapply(A, function(j)
    summary(lm(d$M[, j] ~ d$X + d$Z))$coefficients[2, "t value"], numeric(1))
  expect_equal(tab$t_outcome, t_out, tolerance = 1e-10)
  expect_equal(tab$t_exposure, t_exp, tolerance = 1e-10)
  # The screening p-value is the larger of the two two-sided normal p-values.
  expect_equal(tab$screen_p,
               pmax(2 * pnorm(-abs(t_out)), 2 * pnorm(-abs(t_exp))),
               tolerance = 1e-8)
  # The Bonferroni screen keeps a mediator when its screening p-value is at
  # most error_level / (k * q * log(log(n))), with k = 5 candidates and q = 1.
  keep <- tab$screen_p <= 0.05 / (length(A) * 1 * log(log(nrow(d$X))))
  expect_identical(tab$selected, keep)
  expect_identical(A[keep], as.integer(attr(fit, "selected_mediators")))
  # J_m is sqrt(p) times the summed |t_outcome * t_exposure| of the kept
  # mediators.
  expect_equal(v(fit, "j_pe"),
               sqrt(ncol(d$M)) * sum(abs(t_out * t_exp)[keep]),
               tolerance = 1e-10)
  # tidy() returns the same table as a plain data frame.
  td <- tidy(fit)
  expect_s3_class(td, "data.frame")
  expect_identical(names(td), names(tab))
  for (nm in names(tab)) expect_identical(td[[nm]], tab[[nm]], info = nm)
})

test_that("with two exposures the mediator table keeps the strongest exposure path", quiet_grid({
  d2 <- two_exposure_data("continuous", 0.5)
  # Reversing the sign of the second exposure makes its path t statistics
  # negative, so the strongest path is the largest in magnitude, not the
  # largest signed value.
  X <- d2$X; X[, 2] <- -X[, 2]
  fit <- pe_mediation(X, d2$Y, d2$M, outcome = "continuous",
                      lambda_grid = pinned_grid(200, 30),)
  tab <- pe_mediators(fit)
  A <- tab$mediator
  expect_identical(A, 1:2)
  t_out <- unname(summary(lm(d2$Y ~ d2$M[, A] + X))$coefficients[1 + seq_along(A), "t value"])
  # One row per exposure, one column per candidate mediator.
  t_exp <- vapply(A, function(j)
    summary(lm(d2$M[, j] ~ X))$coefficients[2:3, "t value"], numeric(2))
  expect_true(all(t_exp[2, ] < 0))
  pair_p <- pmax(matrix(2 * pnorm(-abs(t_out)), 2, length(A), byrow = TRUE),
                 2 * pnorm(-abs(t_exp)))
  expect_equal(tab$t_outcome, t_out, tolerance = 1e-10)
  expect_equal(tab$t_exposure, apply(abs(t_exp), 2, max), tolerance = 1e-10)
  expect_equal(tab$screen_p, apply(pair_p, 2, min), tolerance = 1e-6)
  # A mediator is selected when any of its (exposure, mediator) pairs passes
  # the Bonferroni screen over the k * q = 4 pairs; J_m sums the passing
  # pairs only. Here mediator 1 passes through the second exposure alone.
  pass <- pair_p <= 0.05 / (length(A) * 2 * log(log(200)))
  expect_identical(pass, matrix(c(FALSE, TRUE, FALSE, FALSE), 2, 2))
  expect_identical(tab$selected, apply(pass, 2, any))
  expect_equal(v(fit, "j_pe"),
               sqrt(30) * sum(abs(matrix(t_out, 2, 2, byrow = TRUE) * t_exp)[pass]),
               tolerance = 1e-10)
}))

test_that("glance() reports the fit's own statistics, with the PE p-value as pval_pe", quiet_grid({
  fit <- fit_example()
  g <- glance(fit)
  expect_equal(nrow(g), 1L)
  for (t in c("stat_hdmm", "pval_hdmm", "stat_pe", "pval_pe",
              "total_indirect_effect", "n_selected_mediators", "n_observations"))
    expect_identical(g[[t]], v(fit, t), info = t)
  # pval_pe is the upper tail of the PE statistic, not the Wald p-value
  # (about 1e-101 against 0.031 here), so it is compared on the log scale.
  expect_equal(log(g$pval_pe),
               pchisq(g$stat_pe, df = 1, lower.tail = FALSE, log.p = TRUE),
               tolerance = 1e-10)
  expect_equal(g$pval_hdmm, pchisq(g$stat_hdmm, df = 1, lower.tail = FALSE),
               tolerance = 1e-10)
  expect_equal(g$stat_pe, g$stat_hdmm + v(fit, "j_pe"))
  expect_equal(g$stat_pe, 457.768594595678, tolerance = 1e-8)
  expect_equal(g$stat_hdmm, 4.65332089289092, tolerance = 1e-8)

  # With named exposures glance() carries one total indirect effect per
  # exposure, under the exposure's name.
  d2 <- two_exposure_data("continuous", 0.5)
  X <- d2$X; colnames(X) <- c("dose", "age")
  fit2 <- pe_mediation(X, d2$Y, d2$M, outcome = "continuous",
                       lambda_grid = pinned_grid(200, 30),)
  g2 <- glance(fit2)
  for (t in c("total_indirect_effect_dose", "total_indirect_effect_age",
              "stat_pe", "pval_pe"))
    expect_identical(g2[[t]], v(fit2, t), info = t)
  expect_equal(c(g2$total_indirect_effect_dose, g2$total_indirect_effect_age),
               c(0.0304066440499535, 0.218573350987254), tolerance = 1e-8)
}))

# ---------------------------------------------------------------------------
# The BH and BY screens, in cases where they select different mediators.
# ---------------------------------------------------------------------------

test_that("the BH and BY screens reproduce the reference implementation on the continuous example", {
  bh <- fit_example("BH")
  by <- fit_example("BY")
  # At this grid the penalized fit keeps mediators 1 to 5 and the Wald
  # statistic is the same for every method; the screens differ.
  expect_equal(v(bh, "stat_hdmm"), 4.65332089289092, tolerance = 1e-8)
  expect_equal(v(by, "stat_hdmm"), 4.65332089289092, tolerance = 1e-8)
  expect_equal(v(bh, "stat_pe"), 781.826337762419, tolerance = 1e-8)
  expect_equal(v(bh, "j_pe"), 777.173016869529, tolerance = 1e-8)
  expect_equal(log(v(bh, "pval_pe")), log(4.82390829159556e-172), tolerance = 1e-8)
  expect_identical(as.integer(attr(bh, "selected_mediators")), 3:5)
  expect_equal(v(by, "stat_pe"), 457.768594595678, tolerance = 1e-8)
  expect_equal(v(by, "j_pe"), 453.115273702787, tolerance = 1e-8)
  expect_equal(log(v(by, "pval_pe")), log(1.47057154809541e-101), tolerance = 1e-8)
  expect_identical(as.integer(attr(by, "selected_mediators")), 4:5)

  # Both sets follow the step-up rules at error_level / log(log(n)).
  tab <- pe_mediators(by)
  level <- 0.05 / log(log(200))
  expect_identical(tab$mediator[step_up(tab$screen_p, level)], 3:5)
  expect_identical(tab$mediator[step_up(tab$screen_p, level, by = TRUE)], 4:5)
  in_bh <- tab$mediator %in% 3:5
  expect_equal(v(bh, "j_pe"),
               sqrt(100) * sum(abs(tab$t_outcome * tab$t_exposure)[in_bh]),
               tolerance = 1e-10)

  # report_all_methods = TRUE screens with all three methods off the same
  # fit; each pe_selection() row equals the fit run with that method alone.
  all3 <- fit_example("BY", report_all_methods = TRUE)
  sel <- pe_selection(all3)
  expect_identical(all3$value, by$value)
  for (m in c("BH", "BY")) {
    one <- if (m == "BH") bh else by
    row <- sel[sel$method == m, ]
    expect_identical(row$stat_pe, v(one, "stat_pe"), info = m)
    expect_identical(row$j_pe, v(one, "j_pe"), info = m)
    expect_identical(row$pval_pe, v(one, "pval_pe"), info = m)
    expect_identical(row$n_selected, length(attr(one, "selected_mediators")), info = m)
    expect_identical(attr(all3, "selection_by_method")[[m]],
                     attr(one, "selected_mediators"), info = m)
  }
  expect_identical(sel$selected_mediators[sel$method == "BH"], "3, 4, 5")
  expect_identical(sel$selected_mediators[sel$method == "BY"], "4, 5")
})

test_that("BY can select nothing where BH selects, and the PE test then equals the Wald test", quiet_grid({
  d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                               pattern = "contrasting", c1 = 0.75, seed = 113)
  grids <- list(lambda_grid = pinned_grid(200, 60))
  bh <- do.call(pe_mediation, c(list(d$X, d$Y, d$M, outcome = "continuous",
                                     method = "BH"), grids))
  by <- do.call(pe_mediation, c(list(d$X, d$Y, d$M, outcome = "continuous",
                                     method = "BY"), grids))
  # Two candidates survive the penalized fit. The smaller screening p-value
  # (about 0.0142) passes the BH bound, 1/2 of 0.05 / log(log(200)), but not
  # the BY bound, which divides that by 1 + 1/2.
  tab <- pe_mediators(by)
  level <- 0.05 / log(log(200))
  expect_identical(tab$mediator, 1:2)
  expect_identical(tab$mediator[step_up(tab$screen_p, level)], 2L)
  expect_identical(tab$mediator[step_up(tab$screen_p, level, by = TRUE)], integer(0))
  expect_identical(as.integer(attr(bh, "selected_mediators")), 2L)
  expect_identical(as.integer(attr(by, "selected_mediators")), integer(0))
  expect_equal(v(bh, "stat_pe"), 136.405351087556, tolerance = 1e-8)
  expect_equal(v(by, "j_pe"), 0)
  expect_equal(v(by, "stat_pe"), 0.00892114914275966, tolerance = 1e-8)
  expect_identical(v(by, "stat_pe"), v(by, "stat_hdmm"))
  expect_equal(v(by, "pval_pe"), 0.924750241397531, tolerance = 1e-8)
}))

# ---------------------------------------------------------------------------
# Several exposures: the chi-square calibration on q degrees of freedom.
# ---------------------------------------------------------------------------

test_that("with two exposures both p-values are chi-square upper tails on 2 degrees of freedom", quiet_grid({
  # With 2 degrees of freedom the chi-square upper tail is exp(-x / 2), so the
  # p-values are checked against that closed form, on the log scale for the
  # far tail. The moderate p-values also separate df = 2 from df = 1.
  cases <- list(
    list(outcome = "continuous", c1 = 0.3, grid = pinned_grid(200, 30),
         stat_hdmm = 4.10996390822357, stat_pe = 4.10996390822357,
         beta = c(0.0153130955743463, 0.133701896412331), selected = integer(0)),
    list(outcome = "continuous", c1 = 0.5, grid = pinned_grid(200, 30),
         stat_hdmm = 10.8123138682857, stat_pe = 425.555732465452,
         beta = c(0.0304066440499535, 0.218573350987254), selected = 1L),
    list(outcome = "binary", c1 = 0.5, grid = oracle_grid,
         stat_hdmm = 6.19465259798562, stat_pe = 137.519971066375,
         beta = c(-0.0335205644542413, 0.422697761895359), selected = 1L),
    list(outcome = "count", c1 = 0.3, grid = oracle_grid,
         stat_hdmm = 0.202775219981631, stat_pe = 348.498136992233,
         beta = c(-0.0128109087865636, -0.0135999485451871), selected = 4:5))
  fits <- lapply(cases, function(cs) {
    d2 <- two_exposure_data(cs$outcome, cs$c1)
    pe_mediation(d2$X, d2$Y, d2$M, outcome = cs$outcome, lambda_grid = cs$grid)
  })
  for (i in seq_along(cases)) {
    cs <- cases[[i]]; fit <- fits[[i]]
    lab <- paste(cs$outcome, cs$c1)
    expect_identical(v(fit, "df"), 2, info = lab)
    expect_equal(v(fit, "pval_hdmm"), exp(-v(fit, "stat_hdmm") / 2),
                 tolerance = 1e-10, info = lab)
    expect_equal(log(v(fit, "pval_pe")), -v(fit, "stat_pe") / 2,
                 tolerance = 1e-10, info = lab)
    expect_equal(v(fit, "stat_hdmm"), cs$stat_hdmm, tolerance = 1e-8, info = lab)
    expect_equal(v(fit, "stat_pe"), cs$stat_pe, tolerance = 1e-8, info = lab)
    expect_equal(v(fit, "total_indirect_effect_1"), cs$beta[1], tolerance = 1e-8, info = lab)
    expect_equal(v(fit, "total_indirect_effect_2"), cs$beta[2], tolerance = 1e-8, info = lab)
    expect_identical(as.integer(attr(fit, "selected_mediators")), cs$selected, info = lab)
  }
  # At a moderate statistic the reference implementation's own p-value,
  # 1 - pchisq(x, 2), loses nothing to cancellation and agrees.
  expect_equal(v(fits[[1]], "pval_pe"), 0.128095147147381, tolerance = 1e-8)
}))

# ---------------------------------------------------------------------------
# The formula front end with confounders.
# ---------------------------------------------------------------------------

test_that("pe_mediate() with confounders equals pe_mediation() with the confounder matrix", quiet_grid({
  d <- simulate_mediation_data(n = 200, p = 30, outcome = "continuous",
                               pattern = "homogeneous", c1 = 0.5, d = 2,
                               seed = 113)
  df <- data.frame(y = d$Y, x = d$X[, 1], z1 = d$Z[, 1],
                   grp = factor(rep(c("a", "b", "c", "d"), 50)), d$M)
  meds <- grep("^X", names(df), value = TRUE)
  grids <- list(lambda_grid = pinned_grid(200, 30))
  a <- do.call(pe_mediate, c(list(y ~ x, data = df, mediators = meds,
                                  confounders = c("z1", "grp"),
                                  outcome = "continuous"), grids))
  # The factor becomes three indicator columns and no intercept column.
  Z <- model.matrix(~ z1 + grp, data = df)[, -1]
  expect_identical(colnames(Z), c("z1", "grpb", "grpc", "grpd"))
  b <- do.call(pe_mediation, c(list(d$X, d$Y, as.matrix(df[, meds]), Z = Z,
                                    outcome = "continuous"), grids))
  expect_identical(a$term, b$term)
  expect_identical(a$value, b$value)
  expect_identical(attr(a, "selected_mediators"), attr(b, "selected_mediators"))
  expect_identical(tidy(a), tidy(b))
  expect_equal(v(a, "stat_pe"), 230.452405128491, tolerance = 1e-8)
  expect_equal(v(a, "total_indirect_effect"), 0.159536596341202, tolerance = 1e-8)
  expect_identical(as.integer(attr(a, "selected_mediators")), 4L)
  # The confounders matter here: without them mediators 4 and 5 are selected.
  none <- do.call(pe_mediate, c(list(y ~ x, data = df, mediators = meds,
                                     outcome = "continuous"), grids))
  expect_identical(as.integer(attr(none, "selected_mediators")), 4:5)
  # A character confounder is expanded like the factor.
  df_chr <- df; df_chr$grp <- as.character(df_chr$grp)
  chr <- do.call(pe_mediate, c(list(y ~ x, data = df_chr, mediators = meds,
                                    confounders = c("z1", "grp"),
                                    outcome = "continuous"), grids))
  expect_identical(chr$value, a$value)
  # A confounder that is not a column of `data` stops with an error naming it.
  expect_error(pe_mediate(y ~ x, data = df, mediators = meds,
                          confounders = "nope", outcome = "continuous"),
               "nope")
}))
