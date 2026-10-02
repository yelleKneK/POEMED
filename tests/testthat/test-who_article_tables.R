# Cell-by-cell comparison of WHO_mediation_analysis() under its default grids
# with the article's data-analysis tables as printed in the accepted article:
# Table 1 (infant mortality) and Table 2 (prevalence of undernourishment) of
# Yu and Kelley (in press) and Tables S.8 to S.12 of its supplement, which
# give the extended layout (the benchmark p-value, the power-enhanced
# p-value, and the selected set under Bonferroni, BY, and BH) for all five
# outcomes in the global model, each WHO region, each World Bank income
# group, and each region-by-income cell. The printed values are stored in
# who_tables_article_2026-09-30.csv, one row per populated cell, with the
# columns outcome, group, n, article_p_hdmm, article_p_pe, article_selected_pe
# (the Bonferroni set), article_selected_by, and article_selected_bh. In that
# file "<1e-12" is the tables' floor "< 10^-12", "null" is an empty
# selected set, LM and UM are the Lower-middle and Upper-middle income
# groups, and a cell the tables print as "-" (no data) has no row.
#
# For every populated cell the tests check the number of countries, the
# three selected sets as sets, the benchmark p-value to the four decimals
# printed, and the power-enhanced p-value against the printed floor; they
# also check that the cells the package skips at the design stage are
# exactly the ones the tables print as "-", and that the run's messages and
# warnings report those events and nothing else. The known differences
# between the package and the article are listed next, each with its cause,
# and each is asserted in the form the package gives, so a change that
# resolved or moved one shows up here. Each outcome takes 6 to 9 s (about
# 40 s in all), so the fits are skipped on CRAN.

# 1. Cells the package reports as NA (pval_hdmm NA, with a fit-stage note)
#    where the article prints a value. In the first four, two indicators are
#    exactly collinear within the cell (dom_che and ext_che sum to 100 in
#    the Eastern Mediterranean lower-middle-income cell), so the design is
#    rank deficient and the fit fails with a singular-matrix error; the
#    article's values (0.0837, 0.0914, 0.1311, and 0.0389 with
#    pvtd_ncu2021_pc selected) are reproduced by dropping one member of
#    each collinear pair by hand (dom_che in that cell; ext_che, and for
#    LEB also shi_che, in the African upper-middle-income cell), the step
#    ?WHO_mediation_analysis shows. In the last two, the outcome is 2.5 in
#    every country-year (the World Bank's lower reporting bound for the
#    prevalence of undernourishment), which the package refuses because Y
#    is constant; the article prints a p-value of 1 there.
na_cells <- c("imr:EMR-LM", "u5mr:EMR-LM", "leb:AFR-UM", "u5mr:AFR-UM",
              "pou:EUR-High", "pou:WPR-High")

# 2. Benchmark p-values that differ from the printed value in the fourth
#    decimal: the package gives 0.648269 where Table 1 prints 0.6482 (the
#    South-East Asia infant-mortality model) and 0.072579 where Table S.8
#    prints 0.0725 (the high-income Americas cell). Checked to within 1e-4.
fourth_decimal <- c("imr:SEAR", "imr:AMR-High")

# 3. Benchmark p-values the article prints as a rounded power of ten:
#    0.0001, 10^-4, 2 x 10^-8, 10^-6, and 6 x 10^-9, where the package gives
#    3.9e-05, 9.4e-05, 1.3e-08, 2.4e-07, and 5.4e-09 (the power-enhanced
#    p-value is below the printed floor in all five). Checked as below
#    0.00015.
tiny <- c("imr:WPR-High", "u5mr:WPR-High", "leb:EUR", "leb:EUR-High",
          "lbw:WPR-High")

# 4. One selected set: under BH in the high-income under-five mortality
#    model, Table S.11 lists gge_gdp in addition to the package's oops_che,
#    shi_che, vhi_che, and pvtd_ncu2021_pc; the Bonferroni and BY sets
#    agree. Checked as the article's set less that indicator.
bh_extra <- c("u5mr:High" = "gge_gdp")

# 5. Not a difference, but a warning the run emits: in the high-income
#    Eastern Mediterranean undernourishment cell (one country) the
#    penalized fit HBIC selects contains no mediator, so the package reports
#    the p-values of 1 of its empty-fit convention with a warning of class
#    poemed_empty_fit; Table 2 prints 1 as well.
empty_cells <- "pou:EMR-High"

article <- read.csv(test_path("who_tables_article_2026-09-30.csv"),
                    colClasses = "character")
article$n <- as.integer(article$n)
article$cell <- paste(article$outcome, article$group, sep = ":")

# The tables' group labels as WHO_mediation_analysis() spells them.
package_group <- function(g) {
  g <- sub("^LM$", "Lower-middle", g)
  g <- sub("^UM$", "Upper-middle", g)
  g <- sub("-LM$", "-Lower-middle", g)
  sub("-UM$", "-Upper-middle", g)
}

# An selected set as a character vector (empty for the tables' "null" and
# the package's "none"), and as one canonical string, so that two orderings
# of the same set compare equal.
split_set <- function(s) {
  if (is.na(s) || s %in% c("null", "none", "")) return(character(0))
  trimws(strsplit(s, ",")[[1]])
}
as_set <- function(x) {
  vapply(x, function(s) paste(sort(split_set(s)), collapse = ", "), "",
         USE.NAMES = FALSE)
}

# Run one outcome over all groupings in the extended layout, keeping the
# messages and warnings it emits for the checks below.
run_outcome <- function(outcome) {
  msgs <- character()
  warns <- list()
  res <- withCallingHandlers(
    WHO_mediation_analysis(
      outcome, groupings = c("global", "region", "income", "region_income"),
      full_table = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns[[length(warns) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  list(res = res, messages = msgs, warnings = warns)
}

# Compare one outcome's run with its table, cell by cell.
expect_article_table <- function(outcome) {
  a <- article[article$outcome == outcome, ]
  cell <- a$cell
  run <- run_outcome(outcome)
  res <- run$res
  notes <- attr(res, "notes")
  pk <- res[match(package_group(a$group), res$group), ]
  # Two vectors over the table's cells, named by cell so a failure names it.
  pin <- function(pkg, art, i) {
    expect_identical(stats::setNames(pkg[i], cell[i]),
                     stats::setNames(art[i], cell[i]))
  }

  # Every printed cell is a row of the result, and the cells the tables
  # print as "-" are exactly the ones skipped at the design stage.
  expect_false(anyNA(pk$group))
  expect_setequal(res$group[!is.na(res$n_countries)], package_group(a$group))
  expect_setequal(notes$group[notes$stage == "design"],
                  res$group[is.na(res$n_countries)])

  # The number of countries, including in the cells whose fit fails.
  pin(as.integer(pk$n_countries), a$n, seq_along(cell))

  # The NA cells are exactly the documented ones, each from a failed fit,
  # and the empty fit is the documented one.
  is_na <- cell %in% na_cells
  expect_identical(cell[is.na(pk$pval_hdmm)], cell[is_na])
  pe_cols <- c("pval_pe_bonferroni", "pval_pe_bh", "pval_pe_by",
               "selected_bonferroni", "selected_bh", "selected_by")
  expect_true(all(is.na(pk[is_na, pe_cols])))
  expect_setequal(notes$group[notes$stage == "fit"],
                  package_group(a$group[is_na]))
  expect_setequal(notes$group[notes$stage == "empty"],
                  package_group(a$group[cell %in% empty_cells]))
  expect_true(all(notes$stage %in% c("design", "fit", "empty")))

  # The messages and warnings report those events and nothing else: one
  # message naming the skipped cells, one warning naming the failed fits,
  # and one warning of class poemed_empty_fit for the empty fit.
  expect_length(run$messages, 1L)
  expect_match(run$messages, "skipped and reported as NA")
  classes <- vapply(run$warnings, function(w) class(w)[1L], "")
  expect_identical(sum(classes == "poemed_empty_fit"),
                   as.integer(any(cell %in% empty_cells)))
  failed <- run$warnings[classes != "poemed_empty_fit"]
  expect_length(failed, as.integer(any(is_na)))
  for (w in failed) expect_match(conditionMessage(w), "^The fit failed in")

  fitted <- !is_na
  # The selected sets under each method, as sets.
  art_col <- c(bonferroni = "article_selected_pe", by = "article_selected_by",
               bh = "article_selected_bh")
  for (m in names(art_col)) {
    i <- fitted & !(m == "bh" & cell %in% names(bh_extra))
    pin(as_set(pk[[paste0("selected_", m)]]), as_set(a[[art_col[[m]]]]), i)
  }
  for (x in intersect(names(bh_extra), cell)) {
    i <- match(x, cell)
    expect_setequal(split_set(pk$selected_bh[i]),
                    setdiff(split_set(a$article_selected_bh[i]), bh_extra[[x]]))
  }

  # The benchmark p-value to the four decimals printed, except in the
  # documented cells.
  art_p <- suppressWarnings(as.numeric(a$article_p_hdmm))
  std <- fitted & !(cell %in% c(fourth_decimal, tiny))
  pin(sprintf("%.4f", pk$pval_hdmm), sprintf("%.4f", art_p), std)
  for (x in intersect(fourth_decimal, cell)) {
    i <- match(x, cell)
    expect_lt(abs(pk$pval_hdmm[i] - art_p[i]), 1e-4)
  }
  for (x in intersect(tiny, cell)) {
    expect_lt(pk$pval_hdmm[match(x, cell)], 0.00015)
  }

  # The power-enhanced p-value: below the printed floor where the table
  # prints < 10^-12, and equal to the benchmark p-value where no mediator is
  # selected, which is how the tables print it.
  floor <- fitted & a$article_p_pe == "<1e-12"
  expect_identical(cell[floor & !(pk$pval_pe_bonferroni < 1e-12)],
                   character(0))
  pin(pk$pval_pe_bonferroni, pk$pval_hdmm, fitted & !floor)
  # Under BY and BH the tables print it the same way: the floor where the
  # method's set is nonempty, the benchmark p-value where it is empty. The
  # file does not carry those two columns, so the package's values are
  # checked against its own sets.
  for (m in c("by", "bh")) {
    p <- pk[[paste0("pval_pe_", m)]]
    nonempty <- as_set(pk[[paste0("selected_", m)]]) != ""
    expect_identical(cell[fitted & nonempty & !(p < 1e-12)], character(0))
    pin(p, pk$pval_hdmm, fitted & !nonempty)
  }
  invisible(res)
}

test_that("the stored article tables hold every populated cell once", {
  expect_identical(nrow(article), 134L)
  expect_setequal(article$outcome, c("imr", "pou", "leb", "u5mr", "lbw"))
  expect_identical(anyDuplicated(article$cell), 0L)
  expect_true(all(article$n >= 1L))
  expect_true(all(c(na_cells, fourth_decimal, tiny, names(bh_extra),
                    empty_cells) %in% article$cell))
  # A table prints the power-enhanced p-value at the floor exactly where it
  # selects a mediator under Bonferroni.
  expect_identical(article$article_p_pe == "<1e-12",
                   as_set(article$article_selected_pe) != "")
})

test_that("Table 1 and Table S.8: infant mortality", {
  skip_on_cran()
  expect_article_table("imr")
})

test_that("Table 2 and Table S.9: prevalence of undernourishment", {
  skip_on_cran()
  expect_article_table("pou")
})

test_that("Table S.10: life expectancy at birth", {
  skip_on_cran()
  expect_article_table("leb")
})

test_that("Table S.11: under-five mortality", {
  skip_on_cran()
  expect_article_table("u5mr")
})

test_that("Table S.12: low birth weight", {
  skip_on_cran()
  expect_article_table("lbw")
})
