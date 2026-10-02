# Input guards of the fitting functions, the formula interface, and the small
# utilities. Each guard is made to fire here, and each message is checked for
# the name of the argument the user passed, so a bad input stops with a
# message that says what to fix instead of an error from deep inside base R,
# Lapack, or glmnet. Valid calls are checked to return exactly what they
# returned before the guard existed.

guard_data <- function() {
  simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                          pattern = "contrasting", c1 = 1, seed = 113)
}

guard_frame <- function() {
  d <- guard_data()
  set.seed(113)
  df <- data.frame(y = d$Y, x = d$X[, 1], z1 = rnorm(200),
                   g = factor(sample(c("a", "b", "c"), 200, TRUE)), d$M)
  list(df = df, meds = grep("^X", names(df), value = TRUE))
}

# ---- Enumerated arguments name themselves --------------------------------

test_that("an invalid choice names the argument and lists the choices", {
  d <- guard_data()
  expect_error(pe_mediation(d$X, d$Y, d$M, method = "holm"),
               "`method` should be one of \"Bonferroni\", \"BH\", \"BY\"")
  expect_error(pe_mediation(d$X, d$Y, d$M, outcome = "continous"),
               "`outcome` should be one of")
  expect_error(pe_mediation(d$X, d$Y, d$M, outcome = c("continuous", "binary")),
               "`outcome` should be one of .*a single value")
  expect_error(pe_mediation(d$X, d$Y, d$M, outcome = NA), "`outcome`")
  expect_error(pe_mediation_linear(d$X, d$Y, d$M, method = "x"), "`method`")
  expect_error(pe_mediation_logistic(d$X, d$Y, d$M, method = "x"), "`method`")
  expect_error(pe_mediation_poisson(d$X, d$Y, d$M, method = "x"), "`method`")
})

test_that("valid choices behave as under match.arg()", quiet_grid({
  # The default is the first choice, and a unique abbreviation is accepted.
  d <- guard_data()
  expect_identical(pe_mediation(d$X, d$Y, d$M, outcome = "cont", method = "Bonf"),
                   pe_mediation(d$X, d$Y, d$M))
}))

test_that("with several choices allowed, an entry that matches nothing stops", {
  # match.arg() drops such an entry and returns the rest; the shared helper
  # refuses it, and otherwise returns what match.arg() returns.
  ch <- c("alpha", "beta", "betamax")
  for (bad in list(c("alpha", "gamma"), c("alpha", "bet"), c("alpha", NA)))
    expect_error(POEMED:::.poemed_match_arg(bad, ch, "choice", several.ok = TRUE),
                 "`choice` should be one of \"alpha\", \"beta\", \"betamax\", or several of them; you supplied",
                 fixed = TRUE)
  for (ok in list(ch, "alpha", c("betam", "al"), c("beta", "beta"), NULL))
    expect_identical(POEMED:::.poemed_match_arg(ok, ch, "choice", several.ok = TRUE),
                     match.arg(ok, ch, several.ok = TRUE))
})

# ---- Logical flags --------------------------------------------------------

test_that("the logical flags must each be a single TRUE or FALSE", {
  d <- guard_data()
  workers <- list(pe_mediation, pe_mediation_linear, pe_mediation_logistic,
                  pe_mediation_poisson)
  for (w in workers) {
    for (v in list(NA, "yes", c(TRUE, FALSE), 1)) {
      expect_error(w(d$X, d$Y, d$M, report_all_methods = v),
                   "`report_all_methods` must be TRUE or FALSE")
      expect_error(w(d$X, d$Y, d$M, drop_constant = v),
                   "`drop_constant` must be TRUE or FALSE")
      expect_error(w(d$X, d$Y, d$M, scale = v), "`scale` must be TRUE or FALSE")
    }
  }
})

# ---- The outcome ----------------------------------------------------------

test_that("a factor outcome stops instead of being replaced by its level codes", quiet_grid({
  cnt <- example_count
  grid <- seq(0.7, 5, length.out = 15)
  expect_error(pe_mediation(cnt$X, factor(cnt$Y), cnt$M, outcome = "count",
                            lambda_grid = grid),
               "`Y` is a factor")
  expect_error(pe_mediation_poisson(cnt$X, factor(cnt$Y), cnt$M,
                                    lambda_grid = grid),
               "`Y` is a factor")
  d <- guard_data()
  expect_error(pe_mediation_linear(d$X, factor(round(d$Y, 1)), d$M),
               "as.numeric\\(as.character\\(Y\\)\\)")
  expect_error(pe_mediation(example_binary$X, factor(example_binary$Y),
                            example_binary$M, outcome = "binary"),
               "`Y` is a factor")
  # The converted factor gives the numeric fit.
  expect_identical(
    pe_mediation(cnt$X, as.numeric(as.character(factor(cnt$Y))), cnt$M,
                 outcome = "count", lambda_grid = grid)$value,
    pe_mediation(cnt$X, cnt$Y, cnt$M, outcome = "count",
                 lambda_grid = grid)$value)
}))

test_that("a character outcome stops; a logical binary outcome is accepted", quiet_grid({
  d <- guard_data()
  expect_error(pe_mediation(d$X, as.character(d$Y), d$M),
               "`Y` is a character vector")
  b <- example_binary
  expect_identical(pe_mediation(b$X, as.logical(b$Y), b$M, outcome = "binary",
                                lambda_grid = fast_grid),
                   pe_mediation(b$X, b$Y, b$M, outcome = "binary",
                                lambda_grid = fast_grid))
}))

test_that("an atomic outcome that is not numeric stops instead of being coerced", {
  # Complex numbers, dates, and time differences are atomic but not numeric;
  # as.numeric() would drop an imaginary part or read a date as a day count.
  d <- guard_data()
  expect_error(pe_mediation(d$X, complex(real = d$Y), d$M),
               "`Y` must be a numeric vector.", fixed = TRUE)
  expect_error(pe_mediation(d$X, structure(d$Y, class = "Date"), d$M),
               "`Y` must be a numeric vector.", fixed = TRUE)
  expect_error(pe_mediation(d$X, as.difftime(d$Y, units = "secs"), d$M),
               "`Y` must be a numeric vector.", fixed = TRUE)
  expect_error(pe_mediation_linear(d$X, as.difftime(d$Y, units = "secs"), d$M),
               "`Y` must be a numeric vector.", fixed = TRUE)
  cnt <- example_count
  expect_error(pe_mediation(cnt$X, as.difftime(cnt$Y, units = "days"), cnt$M,
                            outcome = "count"),
               "`Y` must be a numeric vector.", fixed = TRUE)
  # A list is refused with the same message.
  expect_error(pe_mediation(d$X, as.list(d$Y), d$M),
               "`Y` must be a numeric vector.", fixed = TRUE)
})

test_that("degenerate X, Y, and M stop with a message naming the argument", {
  d <- guard_data()
  expect_error(pe_mediation(d$X, cbind(d$Y, d$Y), d$M),
               "`Y` must be a single outcome")
  expect_error(pe_mediation(d$X, NULL, d$M), "`Y` must be a non-empty")
  expect_error(pe_mediation(d$X, numeric(0), d$M), "`Y` must be a non-empty")
  expect_error(pe_mediation(NULL, d$Y, d$M), "`X` must be supplied")
  expect_error(pe_mediation(d$X[, 0, drop = FALSE], d$Y, d$M),
               "`X` has no columns")
  expect_error(pe_mediation(d$X, d$Y, NULL), "`M` must be supplied")
  expect_error(pe_mediation(d$X, d$Y, d$M[, 0, drop = FALSE]),
               "`M` must have at least two columns")
})

test_that("a constant outcome is called constant, not binary", {
  d <- guard_data()
  expect_error(pe_mediation_linear(d$X, rep(2.5, 200), d$M), "`Y` is constant")
  expect_error(pe_mediation(d$X, rep(0, 200), d$M), "`Y` is constant")
  expect_error(pe_mediation(d$X, rep(1, 200), d$M, outcome = "binary"),
               "`Y` is constant \\(every value is 1\\)")
  # A two-valued outcome still points to the logistic model.
  expect_error(pe_mediation_linear(d$X, rep(0:1, 100), d$M), "appears to be binary")
})

# ---- Constant mediator columns ---------------------------------------------

test_that("drop_constant = TRUE does not bypass the two-mediator guard", quiet_grid({
  d <- guard_data()
  M1 <- cbind(d$M[, 1], 1, 2)
  expect_error(pe_mediation(d$X, d$Y, M1, drop_constant = TRUE),
               "`M` has only one non-constant column .* at least two")
  expect_error(pe_mediation_logistic(example_binary$X, example_binary$Y,
                                     cbind(example_binary$M[, 1], 0),
                                     drop_constant = TRUE),
               "`M` has only one non-constant column")
  set.seed(113)
  expect_error(pe_mediation_poisson(d$X, rpois(200, 2), M1, drop_constant = TRUE),
               "`M` has only one non-constant column")
  # Two surviving columns still fit, with the drop reported.
  expect_warning(fit <- pe_mediation(d$X, d$Y, cbind(d$M[, 1:2], 1),
                                     drop_constant = TRUE),
                 "dropped: 3")
  expect_s3_class(fit, "poemed_tbl")
}))

# ---- Column names ------------------------------------------------------------

test_that("an NA column name in M falls back to positions without a crash", quiet_grid({
  d <- guard_data()
  plain <- pe_mediation(d$X, d$Y, unname(d$M))
  M_first <- d$M; colnames(M_first) <- c(NA, paste0("m", 2:60))
  fit <- pe_mediation(d$X, d$Y, M_first)
  expect_null(attr(fit, "mediator_names"))
  expect_identical(fit$value, plain$value)
  expect_identical(attr(fit, "selected_mediators"), attr(plain, "selected_mediators"))
  M_last <- d$M; colnames(M_last) <- c(paste0("m", 1:59), NA)
  fit <- pe_mediation(d$X, d$Y, M_last)
  expect_null(attr(fit, "mediator_names"))
  expect_identical(fit$value, plain$value)
}))

# ---- mediation_ar1_cov() ----------------------------------------------------

test_that("mediation_ar1_cov() names p and rho for missing and infinite values", {
  expect_error(mediation_ar1_cov(NA_real_), "`p` must be")
  expect_error(mediation_ar1_cov(NaN), "`p` must be")
  expect_error(mediation_ar1_cov(Inf), "`p` must be")
  expect_error(mediation_ar1_cov(5, NA_real_), "`rho` must be")
  expect_error(mediation_ar1_cov(5, NaN), "`rho` must be")
  expect_error(mediation_ar1_cov(5, -Inf), "`rho` must be")
})

# ---- pe_selection() ------------------------------------------------------------

test_that("pe_selection() stops on an object that is not a fit", quiet_grid({
  d <- guard_data()
  fit <- pe_mediation(d$X, d$Y, d$M)
  not_fits <- list(NULL, 1, "fit", data.frame(a = 1), list(a = 1),
                   pe_mediators(fit))
  for (a in not_fits)
    expect_error(pe_selection(a), "`fit` does not carry a selection record")
  expect_s3_class(pe_selection(fit), "poemed_tbl")
}))

test_that("the pe_selection() example shows screens that differ", quiet_grid({
  set.seed(113)
  d <- simulate_mediation_data(n = 300, p = 100, outcome = "continuous",
                               pattern = "homogeneous", c1 = 1)
  sel <- pe_selection(pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                                   report_all_methods = TRUE))
  expect_identical(sel$method, c("Bonferroni", "BH", "BY"))
  expect_identical(sel$n_selected, c(3L, 5L, 3L))
  expect_identical(sel$selected_mediators, c("3, 4, 5", "2, 3, 4, 5, 14", "3, 4, 5"))
}))

# ---- pe_mediate(): column roles -------------------------------------------------

test_that("pe_mediate() stops when a column is listed in two roles", {
  f <- guard_frame(); df <- f$df; meds <- f$meds
  expect_error(pe_mediate(y ~ x, df, mediators = c("y", meds)),
               "`mediators` must not include the outcome or exposure columns .*: y")
  expect_error(pe_mediate(y ~ x, df, mediators = c("x", meds)),
               "`mediators` must not include .*: x")
  expect_error(pe_mediate(y ~ x, df, mediators = c("z1", meds),
                          confounders = "z1"),
               "`mediators` and `confounders` both list z1")
  expect_error(pe_mediate(y ~ x, df, mediators = c("X1", meds)),
               "`mediators` lists X1 more than once")
  expect_error(pe_mediate(y ~ x, df, mediators = meds, confounders = "x"),
               "`confounders` must not include .*: x")
})

test_that("pe_mediate() stops on a formula or mediator list it cannot use", {
  f <- guard_frame(); df <- f$df; meds <- f$meds
  expect_error(pe_mediate(y ~ 1, df, mediators = meds),
               "`formula` names no exposure")
  expect_error(pe_mediate(y ~ x, df, mediators = character(0)),
               "`mediators` must name at least two")
  expect_error(pe_mediate(y ~ x, df, mediators = "X1"),
               "`mediators` must name at least two")
  expect_error(pe_mediate(y ~ x, df, mediators = 5:10),
               "`mediators` must be a character vector")
  expect_error(pe_mediate(y ~ x, df, mediators = c(meds, NA)),
               "`mediators` must be a character vector")
  expect_error(pe_mediate(y ~ x, df, mediators = meds, confounders = NA_character_),
               "`confounders` must be a character vector")
  expect_error(pe_mediate(y ~ x, df, mediators = meds, confounders = "nope"),
               "confounder column\\(s\\) not in `data`: nope")
  expect_error(pe_mediate(y ~ x, df, mediators = meds, outcome = "cnt"),
               "`outcome` should be one of")
})

test_that("pe_mediate() expands `.` and stops only when it brings in a mediator or confounder", quiet_grid({
  f <- guard_frame(); df <- f$df; meds <- f$meds
  minus <- function(cols)
    stats::as.formula(paste("y ~ . -", paste(cols, collapse = " - ")))
  # A `.` that pulls the mediators in as exposures stops with the role
  # message, which says how to write the formula.
  expect_error(pe_mediate(y ~ ., df, mediators = meds, lambda_grid = fast_grid),
               "`mediators` must not include the outcome or exposure columns named in `formula`: X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, and 50 more. The `.` in `formula` stands for every column",
               fixed = TRUE)
  expect_error(pe_mediate(minus(c("z1", "g", meds[-5])), df, mediators = meds, lambda_grid = fast_grid),
               "`formula`: X5. The `.` in `formula`", fixed = TRUE)
  expect_error(pe_mediate(minus(meds), df, mediators = meds, lambda_grid = fast_grid,
                          confounders = c("z1", "g")),
               "`confounders` must not include .*: z1, g\\. The `\\.` in `formula`")
  # A `.` with the other roles subtracted fits the model the explicit
  # formula fits, as before the guard.
  expect_identical(pe_mediate(minus(c("z1", "g", meds)), df, mediators = meds, lambda_grid = fast_grid),
                   pe_mediate(y ~ x, df, mediators = meds, lambda_grid = fast_grid))
  expect_identical(pe_mediate(minus(meds), df[c("y", "x", meds)],
                              mediators = meds, lambda_grid = fast_grid),
                   pe_mediate(y ~ x, df, mediators = meds, lambda_grid = fast_grid))
  expect_identical(pe_mediate(minus(c("z1", "g", meds)), df, mediators = meds, lambda_grid = fast_grid,
                              confounders = c("z1", "g")),
                   pe_mediate(y ~ x, df, mediators = meds, lambda_grid = fast_grid,
                              confounders = c("z1", "g")))
  # A column only subtracted in `formula` plays no role there.
  expect_identical(pe_mediate(y ~ x - X1, df, mediators = meds, lambda_grid = fast_grid),
                   pe_mediate(y ~ x, df, mediators = meds, lambda_grid = fast_grid))
  # A `.` with every column subtracted leaves no exposure.
  expect_error(pe_mediate(minus(c("x", "z1", "g", meds)), df, mediators = meds, lambda_grid = fast_grid),
               "`formula` names no exposure")
}))

test_that("pe_mediate() treats confounders = character(0) as no confounders", quiet_grid({
  f <- guard_frame()
  expect_identical(
    pe_mediate(y ~ x, f$df, mediators = f$meds, confounders = character(0)),
    pe_mediate(y ~ x, f$df, mediators = f$meds))
}))

test_that("pe_mediate() fits a confounder and matches the matrix interface", quiet_grid({
  f <- guard_frame(); df <- f$df
  fit <- pe_mediate(y ~ x, df, mediators = f$meds, confounders = c("z1", "g"))
  Z <- cbind(z1 = df$z1, gb = as.numeric(df$g == "b"), gc = as.numeric(df$g == "c"))
  ref <- pe_mediation(cbind(x = df$x), df$y, as.matrix(df[f$meds]), Z = Z)
  expect_equal(fit$value, ref$value)
}))

test_that("pe_mediate() stops on a factor or character outcome, naming it", {
  f <- guard_frame(); df <- f$df
  df$y <- factor(round(df$y))
  expect_error(pe_mediate(y ~ x, df, mediators = f$meds),
               "The outcome `y` is a factor")
  df$y <- as.character(f$df$y)
  expect_error(pe_mediate(y ~ x, df, mediators = f$meds),
               "The outcome `y` is a character vector")
})

# ---- pe_mediate(): missing values ------------------------------------------------

test_that("pe_mediate() drops incomplete rows with a message that counts them", quiet_grid({
  f <- guard_frame(); df <- f$df; meds <- f$meds
  complete_fit <- function(rows)
    pe_mediate(y ~ x, df[-rows, ], mediators = meds, confounders = "g",
               lambda_grid = fast_grid)
  # An NA in the outcome: one row dropped, the message says so.
  df_y <- df; df_y$y[3] <- NA
  expect_message(fit <- pe_mediate(y ~ x, df_y, mediators = meds, confounders = "g",
               lambda_grid = fast_grid),
                 "dropped 1 row \\(of 200\\) .*y \\(1\\)")
  expect_identical(fit, complete_fit(3))
  # An NA in a factor confounder no longer surfaces as a `Z` row-count error.
  df_g <- df; df_g$g[5] <- NA
  expect_message(fit <- pe_mediate(y ~ x, df_g, mediators = meds, confounders = "g",
               lambda_grid = fast_grid),
                 "dropped 1 row .*g \\(1\\)")
  expect_identical(fit, complete_fit(5))
  # An NA in a mediator and one in the exposure: both rows, both columns named.
  df_m <- df; df_m$X7[9] <- NA; df_m$x[11] <- NA
  expect_message(fit <- pe_mediate(y ~ x, df_m, mediators = meds, confounders = "g",
               lambda_grid = fast_grid),
                 "dropped 2 rows \\(of 200\\) .*x \\(1\\), X7 \\(1\\)")
  expect_identical(fit, complete_fit(c(9, 11)))
  # Nothing left to analyze.
  df_all <- df; df_all$y <- NA
  expect_error(pe_mediate(y ~ x, df_all, mediators = meds, lambda_grid = fast_grid),
               "no complete cases")
  # Complete data: no message.
  expect_no_message(pe_mediate(y ~ x, df, mediators = meds, lambda_grid = fast_grid))
}))
