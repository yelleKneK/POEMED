v <- function(tab, t) tab$value[tab$term == t]

make_fit <- function(report_all = FALSE) {
  set.seed(113)
  d <- simulate_mediation_data(n = 160, p = 50, outcome = "continuous",
                               pattern = "contrasting", c1 = 1)
  pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
               report_all_methods = report_all)
}

test_that("the fit reports the point estimate of the total indirect effect and no interval", quiet_grid({
  # The article reports no interval for the total indirect effect, and the
  # Wald interval built on the penalized estimate under-covered in finite
  # samples for binary and count outcomes, so the fit carries no interval
  # rows and no confidence-level argument.
  fit <- make_fit()
  expect_identical(fit$term,
                   c("stat_hdmm", "pval_hdmm", "stat_pe", "j_pe", "pval_pe",
                     "total_indirect_effect", "n_selected_mediators", "df",
                     "n_candidate_mediators", "n_observations"))
  expect_true(is.finite(v(fit, "total_indirect_effect")))
  expect_null(attr(fit, "conf_level"))
  set.seed(113)
  d <- simulate_mediation_data(n = 160, p = 50, outcome = "continuous",
                               pattern = "contrasting", c1 = 1)
  expect_error(pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                            conf_level = 0.95), "unused argument")
}))

test_that("pe_mediators returns the per-mediator screen detail", quiet_grid({
  fit <- make_fit()
  tab <- pe_mediators(fit)
  expect_s3_class(tab, "poemed_tbl")
  expect_true(all(c("mediator", "t_outcome", "t_exposure", "screen_p",
                    "selected") %in% names(tab)))
  # The selected rows match the selected set reported on the fit.
  expect_setequal(tab$mediator[tab$selected], attr(fit, "selected_mediators"))
}))

test_that("pe_selection reports all three methods when asked", quiet_grid({
  fit <- make_fit(report_all = TRUE)
  sel <- pe_selection(fit)
  expect_setequal(sel$method, c("Bonferroni", "BH", "BY"))
  expect_true(all(sel$n_selected >= 0))
  # Single-method fit gives a one-row summary.
  expect_equal(nrow(pe_selection(make_fit())), 1L)
}))

test_that("the formula interface matches the matrix interface", quiet_grid({
  set.seed(113)
  d <- simulate_mediation_data(n = 150, p = 40, outcome = "continuous",
                               pattern = "contrasting", c1 = 1)
  df <- data.frame(y = d$Y, x = d$X[, 1], d$M)
  meds <- grep("^X", names(df), value = TRUE)
  a <- pe_mediate(y ~ x, data = df, mediators = meds, outcome = "continuous")
  b <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  expect_equal(v(a, "stat_pe"), v(b, "stat_pe"))
  expect_identical(attr(a, "selected_mediators"), attr(b, "selected_mediators"))
  expect_error(pe_mediate(y ~ x, data = df, mediators = "nope",
                          outcome = "continuous"), "not in `data`")
}))

test_that("plot and broom verbs work on a fit", quiet_grid({
  fit <- make_fit()
  expect_invisible(plot(fit))
  expect_s3_class(generics::tidy(fit), "data.frame")
  g <- generics::glance(fit)
  expect_equal(nrow(g), 1L)
  expect_true("pval_pe" %in% names(g))
}))

test_that("pe_simulation_study stacks patterns", {
  skip_on_cran()
  set.seed(113)
  st <- pe_simulation_study(n = 100, p = 30, outcome = "continuous",
                            c1_grid = c(0, 1), n_rep = 10, seed = 1)
  expect_true("pattern" %in% names(st))
  expect_setequal(unique(st$pattern), c("homogeneous", "contrasting"))
  expect_equal(nrow(st), 4L)
})

test_that("ss_power_pe_mediation recommends a sample size and powers rise with n", {
  skip_on_cran()
  # This small grid stops short of the 0.8 target, so the planner's classed
  # warning is expected and muffled here; test-study_guards.R makes it fire.
  plan <- withCallingHandlers(
    ss_power_pe_mediation(outcome = "continuous", pattern = "contrasting",
                          p = 30, n_grid = c(60, 120), n_rep = 15, seed = 1),
    poemed_target_not_reached = function(w) invokeRestart("muffleWarning"))
  expect_true(all(c("n", "power_pe", "power_hdmm") %in% names(plan)))
  expect_gte(plan$power_pe[plan$n == 120], plan$power_pe[plan$n == 60])
  rn <- attr(plan, "recommended_n")
  expect_true(is.na(rn) || rn %in% plan$n)
  expect_error(ss_power_pe_mediation(outcome = "continuous",
                                     pattern = "contrasting", p = 30,
                                     n_grid = c(60, 120), target_power = 2),
               "in \\(0, 1\\)")
})
