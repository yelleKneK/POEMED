test_that("the display rounds but the stored values keep full precision", quiet_grid({
  set.seed(113)
  d <- simulate_mediation_data(n = 100, p = 40, outcome = "continuous",
                               pattern = "homogeneous", c1 = 0.6)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  # The stored statistic is a full-precision double, not the rounded display.
  s <- fit$value[fit$term == "stat_pe"]
  expect_true(is.numeric(s))
  expect_false(s == signif(s, 1))      # almost surely not a round number
  # The format method returns character columns; the object itself does not.
  disp <- format(fit)
  expect_type(disp$value, "character")
  expect_type(fit$value, "double")
}))

test_that("whole-number rows print without a decimal part", quiet_grid({
  set.seed(113)
  d <- simulate_mediation_data(n = 100, p = 40, outcome = "continuous",
                               pattern = "homogeneous", c1 = 0.6)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  disp <- format(fit)
  expect_identical(disp$value[fit$term == "n_observations"], "100")
  expect_identical(disp$value[fit$term == "df"], "1")
}))

test_that("p-values print to fixed decimals with a floor", {
  # A tiny p-value prints as the floor label, not as 0.0000.
  tab <- data.frame(term = c("pval_pe", "stat_pe"), value = c(1e-15, 12.3))
  tab <- POEMED:::.as_poemed_tbl(tab, p_terms = "pval_pe")
  disp <- format(tab)
  expect_identical(disp$value[1], "< 0.0001")
  # A moderate p-value prints to four decimals.
  tab2 <- data.frame(term = "pval_pe", value = 0.0317)
  tab2 <- POEMED:::.as_poemed_tbl(tab2, p_terms = "pval_pe")
  expect_identical(format(tab2)$value[1], "0.0317")
})

test_that("printing returns the object invisibly and is idempotent in class", quiet_grid({
  set.seed(113)
  d <- simulate_mediation_data(n = 80, p = 30, outcome = "continuous",
                               pattern = "homogeneous")
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  # The printed table is captured and checked rather than written to the log.
  out <- capture.output(res <- withVisible(print(fit)))
  expect_false(res$visible)
  expect_identical(res$value, fit)
  expect_true(any(grepl("Outcome model: continuous", out, fixed = TRUE)))
  twice <- POEMED:::.as_poemed_tbl(fit)
  expect_identical(class(twice), class(fit))   # .as_poemed_tbl is idempotent
}))

# A small fit shared by the display tests below.
display_fit <- function() {
  set.seed(113)
  d <- simulate_mediation_data(n = 80, p = 30, outcome = "continuous",
                               pattern = "homogeneous")
  pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
}

test_that("options(poemed.digits) sets the significant figures of the display", quiet_grid({
  fit <- display_fit()
  stat <- fit$value[fit$term == "stat_pe"]
  default <- format(fit)$value[fit$term == "stat_pe"]
  expect_identical(default, format(signif(stat, 4), digits = 4))
  withr::local_options(poemed.digits = 2L)
  two <- format(fit)$value[fit$term == "stat_pe"]
  expect_identical(two, format(signif(stat, 2), digits = 2))
  expect_false(identical(two, default))
  expect_false(identical(capture.output(print(fit, digits = 4)),
                         capture.output(print(fit))))
  # An explicit argument still overrides the option.
  expect_identical(format(fit, digits = 4)$value[fit$term == "stat_pe"], default)
}))

test_that("digits and digits_p are bounded, with messages that say what to do", quiet_grid({
  fit <- display_fit()
  expect_error(print(fit, digits = 23),
               "`digits` must be a single whole number from 1 to 22", fixed = TRUE)
  expect_error(format(fit, digits = 23), "for example, `digits = 4`", fixed = TRUE)
  expect_error(print(fit, digits_p = 16),
               "`digits_p` must be a single whole number from 1 to 15", fixed = TRUE)
  expect_error(format(fit, digits_p = 0), "`digits_p` must be", fixed = TRUE)
  expect_error(print(fit, digits = NA), "`digits` must be", fixed = TRUE)
  # A bad value set through the option is reported as coming from the option.
  withr::local_options(poemed.digits = 30)
  expect_error(print(fit), "options(poemed.digits = 4)", fixed = TRUE)
  expect_error(format(fit), "getOption(\"poemed.digits\")", fixed = TRUE)
  # The bounds themselves are accepted.
  expect_output(print(fit, digits = 22, digits_p = 15), "Outcome model")
  expect_output(print(fit, digits = 1, digits_p = 1), "Outcome model")
}))

test_that("the p-value floor label follows digits_p", {
  tab <- POEMED:::.as_poemed_tbl(data.frame(term = "pval_pe", value = 1e-6),
                                 p_terms = "pval_pe")
  expect_identical(format(tab)$value, "< 0.0001")
  expect_identical(format(tab, digits_p = 2)$value, "< 0.01")
})

test_that("[ and subset() keep the p-value rows and the footer of the fit", quiet_grid({
  fit <- display_fit()
  cols <- fit[, c("term", "value")]
  expect_s3_class(cols, "poemed_tbl")
  expect_identical(attr(cols, "p_terms"), attr(fit, "p_terms"))
  expect_identical(attr(cols, "tuning"), attr(fit, "tuning"))
  one <- subset(fit, term == "pval_pe")
  expect_identical(format(one)$value, format(fit)$value[fit$term == "pval_pe"])
  out <- capture.output(print(one))
  expect_true(any(grepl("Outcome model: continuous", out, fixed = TRUE)))
  expect_false(any(grepl("e-[0-9]", out)))
  # A single extracted column is still a plain vector.
  expect_identical(fit[, "value"], fit$value)
  expect_identical(fit[fit$term == "df", "value"], 1)
}))

test_that("rbind() keeps a shared footer and drops one the pieces disagree on", quiet_grid({
  fit <- display_fit()
  # Split and put back together: the same table, footer and all.
  back <- rbind(fit[1:5, ], fit[6:10, ])
  expect_equal(back, fit, ignore_attr = "row.names")
  expect_identical(attr(back, "tuning"), attr(fit, "tuning"))
  # Two different fits: no single-fit footer, but the p-value rows keep
  # their format.
  db <- simulate_mediation_data(n = 100, p = 30, outcome = "binary",
                                pattern = "homogeneous", c1 = 1, seed = 113)
  fb <- pe_mediation(db$X, db$Y, db$M, outcome = "binary")
  both <- rbind(fit, fb)
  expect_s3_class(both, "poemed_tbl")
  expect_null(attr(both, "outcome"))
  expect_null(attr(both, "selected_mediators"))
  expect_null(attr(both, "tuning"))
  expect_identical(attr(both, "p_terms"), c("pval_hdmm", "pval_pe"))
  expect_identical(both$value, c(fit$value, fb$value))   # the numbers are untouched
  out <- capture.output(print(both))
  expect_false(any(grepl("Outcome model|Selected mediators|Tuning parameter", out)))
  expect_false(any(grepl("e-[0-9]", out)))
  # A plain data frame carries no footer, so its rows drop the fit's too.
  mixed <- rbind(fit, data.frame(term = "extra", value = 1))
  expect_null(attr(mixed, "outcome"))
  expect_identical(nrow(mixed), nrow(fit) + 1L)
  # rbind.data.frame options pass through.
  expect_identical(rownames(rbind(fit, fb, make.row.names = FALSE)),
                   as.character(seq_len(2L * nrow(fit))))
}))

test_that("pval rows of a long table print as p-values without p_terms", {
  tab <- POEMED:::.as_poemed_tbl(
    data.frame(term = c("stat_pe", "pval_pe", "pval_hdmm"),
               value = c(12.3, 1e-15, 0.03171)))
  expect_null(attr(tab, "p_terms"))
  expect_identical(format(tab)$value, c("12.3", "< 0.0001", "0.0317"))
})
