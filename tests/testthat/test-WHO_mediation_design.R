test_that("the shipped WHO data has the documented structure", {
  data(WHO_health_mediation, package = "POEMED")
  data(WHO_indicator_codebook, package = "POEMED")
  expect_equal(nrow(WHO_health_mediation), 2002L)
  expect_equal(length(unique(WHO_health_mediation$code)), 91L)
  expect_equal(range(WHO_health_mediation$year), c(2000L, 2021L))
  # The 57 indicators are present, complete, and match the codebook.
  ind <- setdiff(names(WHO_health_mediation),
                 c("code", "country", "region", "income", "year",
                   "gdp_growth", "imr", "u5mr", "leb", "lbw", "pou"))
  expect_equal(length(ind), 57L)
  expect_false(anyNA(WHO_health_mediation[, ind]))
  expect_setequal(ind, WHO_indicator_codebook$indicator)
  # Region tallies match the article.
  reg <- tapply(WHO_health_mediation$code, WHO_health_mediation$region,
                function(x) length(unique(x)))
  expect_equal(reg[["AFR"]], 32L); expect_equal(reg[["EUR"]], 9L)
  expect_equal(reg[["SEAR"]], 5L)
})

test_that("WHO_mediation_design assembles a usable design", {
  des <- WHO_mediation_design("imr")
  expect_equal(des$n, 2002L)
  expect_equal(ncol(des$M), 57L)
  expect_equal(ncol(des$X), 1L)
  expect_equal(length(des$Y), des$n)
  # Confounders: 5 region + 3 income indicators + year = 9 columns.
  expect_equal(ncol(des$Z), 9L)
  expect_false(anyNA(des$M))
})

test_that("subsetting drops constant covariates and unused factor levels", {
  eur <- WHO_mediation_design("imr", region = "EUR")
  expect_equal(eur$n, 9L * 22L)
  # Region is constant within EUR, so Z has no region columns; the design is
  # full rank (no all-zero columns from unused factor levels).
  expect_false(any(grepl("region", colnames(eur$Z))))
  expect_equal(qr(cbind(eur$X, eur$Z))$rank, ncol(eur$Z) + 1L)
})

test_that("outcomes with narrower coverage yield fewer observations", {
  expect_equal(WHO_mediation_design("pou")$n, 79L * 21L)   # 2001-2021
  expect_equal(WHO_mediation_design("lbw")$n, 73L * 21L)   # 2000-2020
})

test_that("WHO_mediation_design validates its arguments", {
  expect_error(WHO_mediation_design("imr", region = "XYZ"), "Unknown region")
  expect_error(WHO_mediation_design("imr", income = "rich"), "Unknown income")
  expect_error(WHO_mediation_design("nope"), "`outcome` should be one of",
               fixed = TRUE)
})

test_that("Z always carries the year, and never comes back NULL", {
  # A region-by-income cell holds region and income constant, so only the
  # year remains; a single-year subset would leave it constant, which
  # pe_mediation() rejects (the help page says so).
  cell <- WHO_mediation_design("imr", region = "WPR", income = "High")
  expect_identical(colnames(cell$Z), "year")
  one_year <- WHO_health_mediation[WHO_health_mediation$year == 2010, ]
  des <- WHO_mediation_design("imr", data = one_year)
  expect_true("year" %in% colnames(des$Z))
  expect_error(pe_mediation(des$X, des$Y, des$M, des$Z,
                            outcome = "continuous"),
               "`Z` has constant column")
})

test_that("the documented hand fit of a group with exact partners reproduces the article", quiet_grid({
  skip_on_cran()                       # one penalized fit, 110 x 56
  # ?WHO_mediation_analysis: dom_che and ext_che sum to 100 in this group;
  # dropping dom_che and fitting on the default grid gives the benchmark
  # p-value the article's Table 1 prints for the cell (0.0837; the article's
  # analysis also removed vfa_che, which changes nothing at four decimals),
  # while dropping ext_che instead gives 0.0786.
  cell <- WHO_mediation_design("imr", region = "EMR", income = "Lower-middle")
  expect_equal(range(cell$M[, "dom_che"] + cell$M[, "ext_che"]), c(100, 100))
  keep <- apply(cell$M, 2, sd) > 0 & cell$indicators != "dom_che"
  grid <- seq(0.1, 2, length.out = 100)
  fit <- pe_mediation(scale(cell$X), cell$Y - mean(cell$Y),
                      scale(cell$M[, keep]), Z = cell$Z,
                      outcome = "continuous", scale = FALSE,
                      lambda_grid = grid)
  expect_equal(round(fit$value[fit$term == "pval_hdmm"], 4), 0.0837)
  expect_identical(attr(fit, "selected_mediators"), integer(0))
}))
