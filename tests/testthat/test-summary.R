# Tests for summary.poemed_tbl, the cross-group digest of any grouped POEMED comparison table.

test_that("summary.poemed_tbl surfaces where PE beats the benchmark (default table)", {
  # Fits several models; the release gate and continuous integration run it
  # in full, and CRAN's machines skip it to keep the check short.
  skip_on_cran()
  res <- WHO_mediation_analysis("imr", groupings = c("global", "region"))
  expect_s3_class(res, "poemed_tbl")
  s <- summary(res)
  expect_s3_class(s, "summary.poemed_tbl")
  expect_identical(s$outcome, "imr")
  expect_gte(s$n_pe, s$n_hdmm)                       # PE detects at least as many
  # Under the default grids the global model, SEAR, and WPR are the PE-only
  # detections, as in Table 1 of Yu and Kelley (in press).
  expect_setequal(s$pe_only$group, c("ALL", "SEAR", "WPR"))
  expect_true(all(s$pe_only$pval_hdmm > 0.05))       # benchmark missed them
  expect_true(all(s$pe_only$pval_pe <= 0.05))        # PE caught them
  expect_true("gge_gdp" %in% names(s$top_mediators))
  expect_error(summary(res, alpha_level = 0), "`alpha_level` must be")
  expect_output(print(s), "alpha_level = 0.05")
  # The digest prints p-values as print.poemed_tbl does: four decimals, the
  # floor, and no scientific notation; the object keeps full precision.
  out <- capture.output(print(s))
  expect_true(any(grepl("< 0.0001", out, fixed = TRUE)))
  expect_false(any(grepl("e-[0-9]", out)))
  hdmm <- s$pe_only$pval_hdmm[s$pe_only$group == "SEAR"]
  expect_true(any(grepl(formatC(hdmm, format = "f", digits = 4), out, fixed = TRUE)))
  expect_identical(s$pe_only$pval_pe, res$pval_pe[match(s$pe_only$group, res$group)])
  expect_setequal(names(s), c("outcome", "full_table", "alpha_level", "n_groups",
                              "n_with_data", "n_hdmm", "n_pe", "pe_only",
                              "top_mediators"))
})

# A small grouped comparison table built by hand, so the print method can be
# checked without fitting anything.
fake_comparison <- function(outcome = "imr") {
  tab <- POEMED:::.as_poemed_tbl(data.frame(
    grouping = c("global", "region", "region"),
    group = c("ALL", "AFR", "WPR"), n_countries = c(91L, 32L, 9L),
    pval_hdmm = c(0.0123, 0.4567, 0.872081), pval_pe = c(7.6e-29, 0.4567, 3.8e-83),
    n_selected = c(1L, 0L, 2L),
    selected_mediators = c("m1", "none", "m2, m3"),
    stringsAsFactors = FALSE))
  attr(tab, "outcome") <- outcome
  attr(tab, "full_table") <- FALSE
  tab
}

test_that("the digest header prints alpha_level as supplied and names the outcome", {
  tab <- fake_comparison()
  out <- capture.output(print(summary(tab, alpha_level = 0.9999999)))
  expect_true(any(grepl("alpha_level = 0.9999999", out, fixed = TRUE)))
  expect_identical(out[1L], "POEMED comparison: IMR")
  out <- capture.output(print(summary(tab)))
  expect_true(any(grepl("WPR +9 +0\\.8721 +< 0\\.0001 +m2, m3", out)))
  expect_false(any(grepl("e-[0-9]", out)))
  # Without an outcome label the header does not print a placeholder.
  attr(tab, "outcome") <- NULL
  expect_identical(capture.output(print(summary(tab)))[1L], "POEMED comparison")
  # A column subset keeps the outcome label.
  attr(tab, "outcome") <- "imr"
  sub <- tab[, c("group", "n_countries", "pval_hdmm", "pval_pe", "selected_mediators")]
  expect_identical(summary(sub)$outcome, "imr")
})

test_that("summary of a table that is not a grouped comparison passes ... on", {
  tab <- POEMED:::.as_poemed_tbl(data.frame(c1 = c(0, 0.5, 1),
                                            rejection_pe = c(0.05, 0.6, 1)))
  expect_identical(summary(tab), summary(as.data.frame(tab)))
  expect_identical(summary(tab, digits = 2), summary(as.data.frame(tab), digits = 2))
  expect_false(identical(summary(tab, quantile.type = 1), summary(tab)))
})

test_that("summary.poemed_tbl also works on the extended (full_table) layout", {
  # Fits several models; the release gate and continuous integration run it
  # in full, and CRAN's machines skip it to keep the check short.
  skip_on_cran()
  res <- WHO_mediation_analysis("imr", groupings = c("global", "region"),
                                full_table = TRUE)
  s <- summary(res)
  expect_true(s$full_table)
  expect_setequal(s$pe_only$group, c("ALL", "SEAR", "WPR"))
  # print method runs without error
  expect_output(print(s), "POEMED comparison")
})
