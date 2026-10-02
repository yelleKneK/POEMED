# Guards on the WHO interface. A `data` that is not in the shape of
# WHO_health_mediation, and an enumerated argument that matches nothing,
# stop with a message naming the argument or column. A missing key stops
# only a design that uses its row. The per-group events of
# WHO_mediation_analysis() (a skipped design, indicators dropped for missing
# values, an empty penalized fit) are named to the user and recorded in the
# "notes" attribute, and no reported value changes.

# The table's values without its attributes, for comparing two runs.
who_values <- function(x) {
  attributes(x) <- attributes(x)[c("names", "row.names")]
  class(x) <- "data.frame"
  x
}

test_that("WHO_mediation_design refuses a `data` that is not a data frame", {
  expect_error(WHO_mediation_design("imr", data = NULL),
               "`data` must be a data.frame", fixed = TRUE)
  expect_error(WHO_mediation_design("imr",
                                    data = as.matrix(WHO_health_mediation)),
               "`data` must be a data.frame", fixed = TRUE)
})

test_that("WHO_mediation_design names each required column that is missing", {
  # Before the guard: a missing region or income column returned a design
  # without those confounders, and a missing exposure or outcome gave the
  # false message "Too few complete observations (0)".
  for (col in c("code", "region", "income", "year", "gdp_growth", "imr")) {
    d <- WHO_health_mediation[, names(WHO_health_mediation) != col]
    expect_error(WHO_mediation_design("imr", data = d),
                 paste0("`data` lacks the column `", col, "`"), fixed = TRUE)
  }
})

test_that("WHO_mediation_design needs at least two indicator columns", {
  keep <- c("code", "country", "region", "income", "year", "gdp_growth",
            "imr", "u5mr", "leb", "lbw", "pou", "che_gdp")
  expect_error(WHO_mediation_design("imr", data = WHO_health_mediation[, keep]),
               "fewer than two indicator", fixed = TRUE)
})

test_that("non-numeric model columns stop, named, in both WHO functions", {
  # An extra character column is treated as a mediator; before the guard it
  # made M a character matrix and the analysis returned an NA row.
  d <- WHO_health_mediation
  d$note <- "x"
  expect_error(WHO_mediation_design("imr", data = d),
               "non-numeric column `note`", fixed = TRUE)
  expect_error(WHO_mediation_analysis("imr", groupings = "global", data = d),
               "non-numeric column `note`", fixed = TRUE)
  d <- WHO_health_mediation
  d$year <- as.character(d$year)
  expect_error(WHO_mediation_design("imr", data = d),
               "non-numeric column `year`", fixed = TRUE)
})

test_that("a factor year stops, saying the year is a numeric trend and how to convert it", {
  # A factor year would enter as one indicator per year, a different model.
  d <- WHO_health_mediation
  d$year <- factor(d$year)
  expect_error(WHO_mediation_design("imr", data = d),
               "non-numeric column `year` (factor). The design enters the year as a numeric trend",
               fixed = TRUE)
  expect_error(WHO_mediation_analysis("imr", groupings = "global", data = d),
               "as.numeric(as.character(year))", fixed = TRUE)
  # The conversion the message names gives the design of the numeric year.
  d$year <- as.numeric(as.character(d$year))
  expect_identical(WHO_mediation_design("imr", data = d),
                   WHO_mediation_design("imr"))
})

# The row of Argentina (AMR, Upper-middle income) for 2000.
arg_row <- which(WHO_health_mediation$code == "ARG" &
                   WHO_health_mediation$year == 2000)

test_that("a missing key on a row the design uses stops, naming the country", {
  # Before the guard, a missing region, income group, or year dropped the row
  # from Z alone (Z one row shorter than X), and a missing code was counted
  # as a country.
  d <- WHO_health_mediation
  d$region[arg_row] <- NA
  expect_error(WHO_mediation_design("imr", data = d),
               "missing values in `region` for ARG on rows this design uses;",
               fixed = TRUE)
  d <- WHO_health_mediation
  d$income[arg_row] <- NA
  expect_error(WHO_mediation_design("imr", region = "AMR", data = d),
               "missing values in `income` for ARG on rows this design uses;",
               fixed = TRUE)
  d <- WHO_health_mediation
  d$year[arg_row] <- NA
  expect_error(WHO_mediation_design("imr", income = "Upper-middle", data = d),
               "missing values in `year` for ARG on rows this design uses;",
               fixed = TRUE)
  # With no code the country is named from `country`, and with neither by
  # its row.
  d <- WHO_health_mediation
  d$code[arg_row] <- NA
  expect_error(WHO_mediation_design("imr", region = "AMR", data = d),
               "missing values in `code` for Argentina on rows this design uses;",
               fixed = TRUE)
  d$country <- NULL
  expect_error(WHO_mediation_design("imr", region = "AMR", data = d),
               paste0("for row ", arg_row, " on rows"), fixed = TRUE)
  # Several countries and several keys are listed together.
  d <- WHO_health_mediation
  bra_row <- which(d$code == "BRA")[1L]
  d$region[arg_row] <- NA
  d$year[bra_row] <- NA
  expect_error(WHO_mediation_design("imr", data = d),
               "missing values in `region`, `year` for ARG, BRA on rows this design uses;",
               fixed = TRUE)
})

test_that("a missing key on a row the design does not use leaves the design as it was", {
  # The rows outside the requested region or income group, and the rows
  # dropped for a missing outcome, are not checked, so the design is the
  # one built from the complete data.
  eur <- WHO_mediation_design("imr", region = "EUR")
  d <- WHO_health_mediation
  d$region[arg_row] <- NA
  expect_identical(WHO_mediation_design("imr", region = "EUR", data = d), eur)
  d <- WHO_health_mediation
  d$year[arg_row] <- NA
  expect_identical(WHO_mediation_design("imr", region = "EUR", data = d), eur)
  d <- WHO_health_mediation
  d$income[arg_row] <- NA
  expect_identical(WHO_mediation_design("imr", income = "High", data = d),
                   WHO_mediation_design("imr", income = "High"))
  d <- WHO_health_mediation
  no_lbw <- which(is.na(d$lbw))[1L]
  d$region[no_lbw] <- NA
  d$code[no_lbw] <- NA
  expect_identical(WHO_mediation_design("lbw", data = d),
                   WHO_mediation_design("lbw"))
})

test_that("the analysis skips only the groups whose design uses a row with a missing key", {
  skip_on_cran()                       # fits the six region models twice
  d <- WHO_health_mediation
  d$region[arg_row] <- NA
  expect_message(
    r <- WHO_mediation_analysis("imr", groupings = c("global", "region"),
                                data = d),
    "1 group skipped and reported as NA: ALL (`data` has missing values in `region` for ARG on rows this design uses;",
    fixed = TRUE)
  nt <- attr(r, "notes")
  expect_identical(nt$group, "ALL")
  expect_identical(nt$stage, "design")
  expect_true(is.na(r$pval_pe[r$group == "ALL"]))
  # The row sits in no region, so every region is fitted as on the data
  # without it.
  regions <- r[r$grouping == "region", ]
  rownames(regions) <- NULL
  expect_identical(
    who_values(regions),
    who_values(WHO_mediation_analysis("imr", groupings = "region",
                                      data = WHO_health_mediation[-arg_row, ])))
})

test_that("a character income column is accepted", {
  # read.csv() and dplyr give a character column, which has no levels; the
  # income check used to call every valid group "Unknown".
  d <- WHO_health_mediation
  d$income <- as.character(d$income)
  a <- WHO_mediation_design("imr", income = "Low", data = d)
  b <- WHO_mediation_design("imr", income = "Low")
  expect_identical(a$n_countries, 16L)
  expect_identical(a[c("X", "Y", "M", "Z", "n")], b[c("X", "Y", "M", "Z", "n")])
  expect_identical(WHO_mediation_design("imr", data = d)$Z,
                   WHO_mediation_design("imr")$Z)
  expect_error(WHO_mediation_design("imr", income = "rich", data = d),
               "Unknown income group(s): rich.", fixed = TRUE)
  # A region that lacks an income group still reports the empty subset.
  expect_error(WHO_mediation_design("imr", region = "AFR", income = "High"),
               "Too few complete observations (0)", fixed = TRUE)
})

test_that("the income grouping runs on a character income column", {
  skip_on_cran()                       # fits four penalized models
  d <- WHO_health_mediation
  d$income <- as.character(d$income)
  expect_identical(
    who_values(WHO_mediation_analysis("imr", groupings = "income", data = d)),
    who_values(WHO_mediation_analysis("imr", groupings = "income")))
})

test_that("enumerated arguments name themselves when they match nothing", {
  expect_error(WHO_mediation_design("xyz"),
               "`outcome` should be one of \"imr\", \"u5mr\"", fixed = TRUE)
  expect_error(WHO_mediation_design(c("imr", "leb")), "`outcome`",
               fixed = TRUE)
  expect_error(WHO_mediation_design(NA), "`outcome`", fixed = TRUE)
  expect_error(WHO_mediation_design(1), "`outcome`", fixed = TRUE)
  expect_error(WHO_mediation_analysis("xyz"), "`outcome` should be one of",
               fixed = TRUE)
  expect_error(WHO_mediation_analysis("imr", groupings = "foo"),
               "`groupings` should be one of", fixed = TRUE)
  expect_error(WHO_mediation_analysis("imr", groupings = character(0)),
               "`groupings`", fixed = TRUE)
  # match.arg() dropped an entry that matched nothing and ran the rest.
  expect_error(WHO_mediation_analysis("imr", groupings = c("global", "regoin")),
               "you supplied c(\"global\", \"regoin\")", fixed = TRUE)
})

test_that("enumerated arguments accept what match.arg() accepted", {
  expect_identical(WHO_mediation_design()$outcome, "imr")
  expect_identical(WHO_mediation_design("u5")$outcome, "u5mr")
  groupings <- c("global", "region", "income", "region_income")
  expect_identical(
    POEMED:::.WHO_match_arg(c("global", "inc"), groupings, "groupings",
                            several.ok = TRUE),
    c("global", "income"))
  # An ambiguous abbreviation ("reg": region or region_income) matches
  # nothing, as in match.arg(), and now stops.
  expect_error(
    POEMED:::.WHO_match_arg(c("global", "reg"), groupings, "groupings",
                            several.ok = TRUE),
    "`groupings`", fixed = TRUE)
  expect_identical(POEMED:::.WHO_match_arg(NULL, groupings, "groupings"),
                   "global")
  # A missing entry is refused too, where match.arg() would drop it.
  expect_error(
    POEMED:::.WHO_match_arg(c("global", NA), groupings, "groupings",
                            several.ok = TRUE),
    "`groupings` should be one of", fixed = TRUE)
  expect_identical(
    POEMED:::.WHO_match_arg(groupings, groupings, "groupings",
                            several.ok = TRUE),
    groupings)
})

test_that("indicators dropped for missing values are named and recorded", {
  skip_on_cran()                       # fits the global model three times
  d <- WHO_health_mediation
  d$che_gdp[5] <- NA                   # Argentina, 2004
  expect_message(
    r <- WHO_mediation_analysis("imr", groupings = "global", data = d),
    "dropped before fitting in 1 group: ALL (che_gdp)", fixed = TRUE)
  nt <- attr(r, "notes")
  expect_identical(nt$group, "ALL")
  expect_identical(nt$stage, "missing")
  expect_match(nt$message, "che_gdp", fixed = TRUE)
  # The same indicator is dropped as before, so the row equals a fit
  # without che_gdp.
  without <- WHO_health_mediation[, names(WHO_health_mediation) != "che_gdp"]
  expect_identical(
    who_values(r),
    who_values(WHO_mediation_analysis("imr", groupings = "global",
                                      data = without)))
  # An infinite value is treated the same way.
  d$che_gdp[5] <- Inf
  expect_message(WHO_mediation_analysis("imr", groupings = "global", data = d),
                 "ALL (che_gdp)", fixed = TRUE)
})

test_that("a missing indicator value on a row the design drops is not reported", {
  skip_on_cran()
  d <- WHO_health_mediation
  d$che_gdp[which(is.na(d$lbw))[1L]] <- NA
  expect_no_message(r <- WHO_mediation_analysis("lbw", groupings = "global",
                                                data = d))
  expect_null(attr(r, "notes"))
  expect_identical(who_values(r),
                   who_values(WHO_mediation_analysis("lbw",
                                                     groupings = "global")))
})

test_that("an empty penalized fit is named and recorded at any core count", {
  skip_on_cran()
  # The undernourishment data for the one Eastern Mediterranean high-income
  # country (21 years): the income grouping skips the three empty income
  # groups and fits "High", whose penalized selection is empty.
  emr_high <- WHO_health_mediation[WHO_health_mediation$region == "EMR" &
                                     WHO_health_mediation$income == "High", ]
  runs <- lapply(c(1L, 2L), function(k) {
    w <- suppressMessages(expect_warning(
      r <- WHO_mediation_analysis("pou", groupings = "income",
                                  data = emr_high, cores = k),
      "The penalized fit that HBIC selected contains no mediator in 1 group: High (lambda = 2).",
      fixed = TRUE, class = "poemed_empty_fit"))
    # The single-fit warning's wording: the reporting convention is the
    # package's, and there is no advice to extend the grid.
    expect_match(conditionMessage(w), "POEMED's reporting convention", fixed = TRUE)
    expect_no_match(conditionMessage(w), "smaller values", fixed = TRUE)
    nt <- attr(r, "notes")
    expect_identical(nt$group[nt$stage == "empty"], "High")
    expect_match(nt$message[nt$stage == "empty"],
                 "The penalized fit that HBIC selected (lambda = 2) contains no mediator",
                 fixed = TRUE)
    # The row keeps the convention values, which the article prints too.
    high <- r[r$group == "High", ]
    expect_equal(c(high$pval_hdmm, high$pval_pe), c(1, 1))
    expect_identical(high$selected_mediators, "none")
    r
  })
  expect_identical(runs[[1L]], runs[[2L]])
})

test_that("the forked run names the empty undernourishment cell EMR-High", {
  skip_on_cran()                       # the full region-by-income table
  # With cores > 1 the fits run in forked children, where the single-fit
  # warning used to be lost with nothing to mark the row.
  suppressWarnings(suppressMessages(expect_warning(
    r <- WHO_mediation_analysis("pou", groupings = "region_income",
                                cores = 2L),
    "EMR-High", class = "poemed_empty_fit")))
  nt <- attr(r, "notes")
  expect_true("EMR-High" %in% nt$group[nt$stage == "empty"])
  expect_equal(r$pval_pe[r$group == "EMR-High"], 1)
})

test_that("a missing region or income that does not enter Z leaves the design unchanged", {
  # Region and income enter the confounders only where they vary among the
  # rows a design keeps. Within one region, a missing region on a kept row
  # changes nothing; where regions vary, the same gap stops the design.
  w <- WHO_health_mediation[WHO_health_mediation$region == "EUR", ]
  d0 <- WHO_mediation_design("imr", income = "High", data = w)
  kept <- which(w$income == "High" & !is.na(w$imr) & !is.na(w$gdp_growth))
  w2 <- w
  w2$region[kept[1]] <- NA
  d1 <- WHO_mediation_design("imr", income = "High", data = w2)
  expect_identical(d1[c("X", "Y", "M", "Z")], d0[c("X", "Y", "M", "Z")])
  full <- WHO_health_mediation
  hit <- which(full$income == "High" & !is.na(full$imr) & !is.na(full$gdp_growth))[1]
  full$region[hit] <- NA
  expect_error(WHO_mediation_design("imr", income = "High", data = full),
               "missing values in `region`")
})
