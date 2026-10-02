# Benchmark numerical-correctness tests. WHO_health_mediation is shipped as a
# benchmark data set and WHO_mediation_analysis() as the benchmark fit. Under
# the function's default tuning grids these tests pin the fit to Table 1
# (infant mortality) and Table 2 (prevalence of undernourishment) of Yu and
# Kelley (in press) as printed in the accepted article: the number of
# countries, the benchmark p-value to the four decimals the tables print, the
# power-enhanced p-value below the printed floor (< 10^-12) where a mediator
# is selected and equal to the benchmark p-value where none is, and the
# complete selected set in every row (a spurious extra mediator fails here). A
# change to the estimation internals that silently moved the published
# findings would fail here. The region-by-income rows of both tables and
# Tables S.8 to S.12 of the supplement are checked cell by cell in
# test-who_article_tables.R.

# The tables print benchmark p-values to four decimals.
printed <- function(p) sprintf("%.4f", p)

# Pin one row of a table. `p_hdmm` is the benchmark p-value as the table
# prints it, or a number where the package's own value is pinned instead;
# `selected` is the printed selected set (empty where the table prints none).
expect_row <- function(res, group, n, p_hdmm, selected = character(0)) {
  r <- res[res$group == group, ]
  expect_equal(nrow(r), 1L, label = paste(group, "rows"))
  expect_equal(r$n_countries, n, label = paste(group, "n_countries"))
  if (is.character(p_hdmm)) {
    expect_identical(printed(r$pval_hdmm), p_hdmm,
                     label = paste(group, "benchmark p-value as printed"))
  } else {
    expect_equal(r$pval_hdmm, p_hdmm, tolerance = 1e-7,
                 label = paste(group, "benchmark p-value"))
  }
  if (length(selected)) {
    # The table prints the power-enhanced p-value as < 10^-12 ...
    expect_lt(r$pval_pe, 1e-12, label = paste(group, "PE p-value"))
    expect_setequal(strsplit(r$selected_mediators, ", ")[[1]], selected)
  } else {
    # ... and, where no mediator is selected, as the benchmark p-value.
    expect_identical(r$pval_pe, r$pval_hdmm,
                     label = paste(group, "PE p-value"))
    expect_identical(r$n_selected, 0L, label = paste(group, "n_selected"))
    expect_identical(r$selected_mediators, "none",
                     label = paste(group, "selected_mediators"))
  }
}

test_that("Table 1 of the article: infant mortality under the default grid", {
  skip_on_cran()                     # 11 penalized fits, about 4 s
  res <- WHO_mediation_analysis("imr",
                                groupings = c("global", "region", "income"))
  expect_identical(res$group,
                   c("ALL", "AFR", "AMR", "EMR", "EUR", "SEAR", "WPR",
                     "Low", "Lower-middle", "Upper-middle", "High"))
  # Global model: general government expenditure as a percent of GDP.
  expect_row(res, "ALL", 91L, "0.0644", "gge_gdp")
  # WHO regions. Africa, the Americas, the Eastern Mediterranean, and Europe:
  # no mediator selected (in Europe the benchmark test rejects on its own;
  # its p-value, 9.4e-05, prints as 0.0001).
  expect_row(res, "AFR", 32L, "0.4434")
  expect_row(res, "AMR", 28L, "0.4974")
  expect_row(res, "EMR", 8L, "0.0261")
  expect_row(res, "EUR", 9L, "0.0001")
  # South-East Asia: external health expenditure in constant 2021 US dollars
  # per capita. The package's benchmark p-value is 0.648269, which Table 1
  # prints as 0.6482 (the one rounding difference among the rows pinned
  # here); the row is pinned to the package's own value.
  expect_row(res, "SEAR", 5L, 0.64826907, "ext_usd2021_pc")
  # Western Pacific: the benchmark Wald test sees nothing, the PE test
  # detects compulsory health insurance as a percent of current health
  # expenditure (the article's headline pattern).
  expect_gt(res$pval_hdmm[res$group == "WPR"], 0.5)
  expect_row(res, "WPR", 9L, "0.8681", "chi_che")
  # World Bank income groups. Low: domestic private health expenditure in
  # million constant 2021 US dollars; High: out-of-pocket spending and
  # social health insurance as percents of current health expenditure; the
  # two middle groups: none.
  expect_row(res, "Low", 16L, "0.3091", "pvtd_usd2021")
  expect_row(res, "Lower-middle", 27L, "0.9543")
  expect_row(res, "Upper-middle", 26L, "0.4481")
  expect_row(res, "High", 22L, "0.0389", c("oops_che", "shi_che"))
})

test_that("Table 2 of the article: the global undernourishment row", {
  skip_on_cran()                     # one fit, under a second, but the
                                     # printed pins stay off CRAN's platforms
  res <- WHO_mediation_analysis("pou", groupings = "global")
  # Social health insurance as a percent of current health expenditure: the
  # benchmark Wald test misses it, the PE test detects it.
  expect_row(res, "ALL", 79L, "0.1626", "shi_che")
})

test_that("Table 2 of the article: the region and income rows", {
  skip_on_cran()                     # 10 penalized fits, about 3 s
  res <- WHO_mediation_analysis("pou", groupings = c("region", "income"))
  expect_identical(res$group,
                   c("AFR", "AMR", "EMR", "EUR", "SEAR", "WPR",
                     "Low", "Lower-middle", "Upper-middle", "High"))
  # Africa: external health expenditure in million constant 2021 US dollars,
  # which the benchmark test misses; Europe: general government expenditure
  # as a percent of GDP.
  expect_row(res, "AFR", 32L, "0.9361", "ext_usd2021")
  expect_row(res, "AMR", 23L, "0.3526")
  expect_row(res, "EMR", 5L, "0.4431")
  expect_row(res, "EUR", 9L, "0.0778", "gge_gdp")
  expect_row(res, "SEAR", 4L, "0.2164")
  expect_row(res, "WPR", 6L, "0.8266")
  # Low income: domestic private health expenditure in million constant 2021
  # US dollars; high income: current health expenditure as a percent of GDP
  # (its benchmark p-value, 7.6e-05, prints as 0.0001).
  expect_row(res, "Low", 16L, "0.0133", "pvtd_usd2021")
  expect_row(res, "Lower-middle", 26L, "0.6506")
  expect_row(res, "Upper-middle", 20L, "0.8245")
  expect_row(res, "High", 17L, "0.0001", "che_gdp")
})
