# The non-key, non-outcome columns of WHO_health_mediation are exactly the 57
# candidate mediators. Computed once here so the design builder and the data
# documentation agree on the set.
#' @keywords internal
#' @noRd
.WHO_indicator_names <- function(data) {
  setdiff(names(data),
          c("code", "country", "region", "income", "year",
            "gdp_growth", "imr", "u5mr", "leb", "lbw", "pou"))
}

# The enumerated arguments of the WHO functions (`outcome`, `groupings`),
# resolved by the package's shared .poemed_match_arg(): valid input behaves
# as under match.arg(), and anything else stops with a message that names
# the argument. With several.ok = TRUE an entry that matches nothing stops
# too, where match.arg() would drop it silently (groupings = c("global",
# "regoin") would run the global model alone). Not exported.
#' @keywords internal
#' @noRd
.WHO_match_arg <- function(arg, choices, name, several.ok = FALSE) {
  .poemed_match_arg(arg, choices, name, several.ok = several.ok)
}

# The data must be a data.frame in the shape of WHO_health_mediation: the key
# and confounder columns, the exposure, the requested outcome, and at least
# two indicator columns, with every column that enters a matrix numeric.
# Checked by WHO_mediation_design() (the engine) and once up front by
# WHO_mediation_analysis(), so a wrong `data` stops with a message naming
# the column instead of surfacing as "subscript out of bounds", a false
# "Too few complete observations (0)", a character mediator matrix, or a
# design that silently lacks the region or income indicators. Region and
# income may be character or factor. Missing keys are not checked here:
# they matter only on the rows a design keeps (see .WHO_check_keys()).
# Not exported.
#' @keywords internal
#' @noRd
.WHO_validate_data <- function(data, outcome) {
  if (!is.data.frame(data))
    stop("`data` must be a data.frame in the shape of WHO_health_mediation; ",
         "pass WHO_health_mediation or a data frame with the same columns.",
         call. = FALSE)
  ticks <- function(x) paste0("`", x, "`", collapse = ", ")
  need <- c("code", "region", "income", "year", "gdp_growth", outcome)
  miss <- setdiff(need, names(data))
  if (length(miss))
    stop("`data` lacks the column", if (length(miss) > 1L) "s" else "", " ",
         ticks(miss), " that WHO_health_mediation carries; ",
         if (length(miss) > 1L)
           "add them, or give the matching columns those names." else
           "add it, or give the matching column that name.", call. = FALSE)
  indicators <- .WHO_indicator_names(data)
  if (length(indicators) < 2L)
    stop("`data` carries fewer than two indicator (mediator) columns; ",
         "supply at least two (every column other than the key, exposure, ",
         "and outcome columns of WHO_health_mediation is treated as a ",
         "candidate mediator).", call. = FALSE)
  # The year enters the confounders as a numeric trend. A factor or
  # character year would instead enter as one indicator per year, a
  # different model, so it stops with the conversion to use: as.numeric()
  # alone would turn a factor into its level codes, not its years.
  if (!is.numeric(data$year))
    stop("`data` has non-numeric column `year` (", class(data$year)[1L],
         "). The design enters the year as a numeric trend, not as one ",
         "indicator per year, so `year` must hold the calendar year as a ",
         "number; convert a factor or character column with ",
         "as.numeric(as.character(year)), which reads the labels rather ",
         "than a factor's level codes.", call. = FALSE)
  numeric_cols <- c("gdp_growth", outcome, indicators)
  bad <- numeric_cols[!vapply(data[numeric_cols], is.numeric, logical(1))]
  if (length(bad))
    stop("`data` has non-numeric column", if (length(bad) > 1L) "s" else "",
         " ", ticks(bad), ". The exposure, the outcome, and every ",
         "indicator column must be numeric (every column other than the key, ",
         "exposure, and outcome columns is treated as a candidate mediator); ",
         "convert ", if (length(bad) > 1L) "them" else "it",
         " with as.numeric() or remove ",
         if (length(bad) > 1L) "them" else "it", ".", call. = FALSE)
  invisible(TRUE)
}

# Every row a design keeps needs its country code, region, income group, and
# year. A missing region, income group, or year would drop the row from the
# confounder matrix alone (model.matrix() omits it), leaving Z shorter than
# X, Y, and M, and a missing code would be counted as a country. Checked by
# WHO_mediation_design() on the rows it keeps, after the region and income
# subset and the drop of rows with no outcome or exposure, so a gap in a row
# the design does not use leaves the design as it was; the message names
# the countries (by `code`, else by `country`, else by row). Not exported.
#' @keywords internal
#' @noRd
.WHO_check_keys <- function(df) {
  # Code and year are needed on every kept row. Region and income matter
  # only where they vary among the kept rows, because only then do they
  # enter the confounders Z (a constant factor is left out of Z); a gap in a
  # key that does not enter Z leaves the design exactly as it was.
  varies <- function(v) length(unique(stats::na.omit(as.character(v)))) > 1L
  keys <- c("code", if (varies(df$region)) "region",
            if (varies(df$income)) "income", "year")
  na <- is.na(df[keys])
  if (!any(na)) return(invisible(TRUE))
  rows <- which(rowSums(na) > 0L)
  who <- as.character(df$code[rows])
  if ("country" %in% names(df)) {
    fill <- is.na(who)
    who[fill] <- as.character(df$country[rows][fill])
  }
  fill <- is.na(who)
  who[fill] <- paste("row", rownames(df)[rows][fill])
  who <- unique(who)
  shown <- utils::head(who, 10L)
  cols <- keys[colSums(na) > 0L]
  stop("`data` has missing values in ", paste0("`", cols, "`", collapse = ", "),
       " for ", paste(shown, collapse = ", "),
       if (length(who) > 10L) sprintf(", and %d more", length(who) - 10L),
       " on rows this design uses; every row the design uses needs its ",
       "country code, region, income group, and year. Remove those rows or ",
       "fill them in.", call. = FALSE)
}

#' Build a Mediation Design From the WHO Health-Expenditure Data
#'
#' Assembles the exposure, outcome, mediator, and confounder matrices for
#' one of the analyses in Yu and Kelley (in press) from
#' [WHO_health_mediation], ready to pass to [pe_mediation()]. The exposure
#' is GDP per capita growth, the mediators are the 57 health-expenditure
#' indicators, and the outcome is the chosen health indicator. Optionally
#' restrict to a WHO region and/or a World Bank income group, as in the
#' article's region-specific, income-specific, and region-by-income
#' models.
#'
#' @details
#' Rows with a missing outcome or missing exposure are dropped (the
#' article handles each outcome's coverage separately, which is why the
#' low-birthweight and undernourishment outcomes have fewer observations).
#' The confounder matrix \eqn{Z}{Z} contains the WHO region and income group
#' (as indicator variables) and the year (numeric), the covariate set of
#' the article. The region and income indicators enter only while they
#' still vary after subsetting: a region-specific model drops the (now
#' constant) region, and unused factor levels are dropped so no all-zero
#' indicator columns remain. The year always enters, so a subset must span
#' at least two years; a single-year subset leaves a constant year column,
#' which [pe_mediation()] rejects.
#'
#' The indicator columns are passed through as they are. A missing value
#' in one of them makes [pe_mediation()] stop, since the method has no
#' missing-data mechanism; [WHO_mediation_analysis()] instead drops such an
#' indicator within the affected group and names it in a message.
#'
#' To fit with the preprocessing of the analysis behind the article's
#' tables, standardize the exposure and mediators and center the outcome
#' while leaving the confounders (the region and income indicators and the
#' year) unscaled. The article's text says that covariates are
#' standardized; the unscaled confounders are a choice of that analysis.
#' The grid below, `seq(0.1, 2, length.out = 100)`, is the grid behind the
#' article's Tables 1, 2, and S.8 to S.12 and the default of
#' [WHO_mediation_analysis()]; the call gives the first row of Table 1
#' (benchmark p-value 0.0644, `gge_gdp` selected):
#' \preformatted{
#'   des <- WHO_mediation_design("imr")
#'   pe_mediation(scale(des$X), des$Y - mean(des$Y), scale(des$M),
#'                Z = des$Z, outcome = "continuous", scale = FALSE,
#'                lambda_grid = seq(0.1, 2, length.out = 100))
#' }
#' HBIC chooses 0.1, the smallest value of that grid, for this fit, as it
#' does in about half of the article's cells, so the fit warns (class
#' `poemed_grid_boundary`); [WHO_mediation_analysis()] keeps the article's
#' grid and does not raise the warning (see its Details).
#'
#' The simpler call `pe_mediation(des$X, des$Y, des$M, des$Z, outcome =
#' "continuous")` standardizes \eqn{X}{X}, \eqn{M}{M}, and \eqn{Z}{Z}
#' and centers \eqn{Y}{Y}, and so fits a different model: scaling the
#' confounder columns changes both the benchmark Wald test and the penalized
#' selection, so its p-values and its selected set can differ from those of
#' the call above. [WHO_mediation_analysis()] applies the preprocessing and
#' the grid of the article's analysis for you.
#'
#' @param outcome Which health outcome to use as \eqn{Y}{Y}: `"imr"` (infant
#'   mortality rate), `"u5mr"` (under-five mortality rate), `"leb"` (life
#'   expectancy at birth), `"lbw"` (prevalence of low birthweight), or
#'   `"pou"` (prevalence of undernourishment).
#' @param region Optional character vector of WHO regions to keep
#'   (`"AFR"`, `"AMR"`, `"EMR"`, `"EUR"`, `"SEAR"`, `"WPR"`). Default
#'   `NULL` keeps all regions.
#' @param income Optional character vector of income groups to keep
#'   (`"Low"`, `"Lower-middle"`, `"Upper-middle"`, `"High"`). Default
#'   `NULL` keeps all groups.
#' @param data The source data frame, in the shape of
#'   [WHO_health_mediation]: the columns `code`, `region`, `income`, `year`,
#'   `gdp_growth`, and the chosen outcome, plus at least two indicator
#'   columns. Every column other than the key, exposure, and outcome
#'   columns of [WHO_health_mediation] is treated as a candidate mediator,
#'   so the year, the exposure, the outcome, and every indicator column
#'   must be numeric; `region` and `income` may be character or factor.
#'   Every row the design uses needs its country `code` and `year`, and
#'   its `region` and `income` wherever those vary among the rows kept (that
#'   is, wherever they enter the confounders); a row missing one of them
#'   stops the design with a message naming the country. Rows outside the
#'   requested region or income group, or dropped for a missing outcome or
#'   exposure, are not checked.
#'   A `data` that breaks these rules stops with a message naming the
#'   column. Defaults to [WHO_health_mediation].
#'
#' @return A list with `X` (exposure, \eqn{n \times 1}{n x 1}), `Y` (outcome
#'   vector), `M` (the \eqn{n \times 57}{n x 57} mediator matrix), `Z` (the
#'   confounder matrix: the region and income indicators that still vary
#'   after subsetting, and the year), and the metadata `outcome`, `n`
#'   (country-year observations), `n_countries` (unique countries),
#'   `indicators`, `region`, and `income`.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @references
#' Yu, X., & Kelley, K. (in press). Power Enhancement in
#' High-Dimensional Heterogeneous Mediation Analysis. \emph{Journal of the
#' American Statistical Association}.
#'
#' @seealso [pe_mediation()], [WHO_health_mediation],
#'   [WHO_indicator_codebook].
#'
#' @family WHO mediation data
#'
#' @examples
#' # Global infant-mortality model: GDP growth -> health spending -> IMR,
#' # with the preprocessing and the grid of the article's Table 1; the fit
#' # selects gge_gdp. HBIC chooses 0.1, the smallest value of that grid, so
#' # the fit warns (class poemed_grid_boundary); WHO_mediation_analysis()
#' # keeps the article's grid and does not raise the warning.
#' des <- WHO_mediation_design("imr")
#' c(n = des$n, mediators = ncol(des$M))
#' fit <- pe_mediation(scale(des$X), des$Y - mean(des$Y), scale(des$M),
#'   Z = des$Z, outcome = "continuous", scale = FALSE,
#'   lambda_grid = seq(0.1, 2, length.out = 100))
#' fit
#' # The selected mediator's full name:
#' WHO_indicator_codebook$description[
#'   match(des$indicators[attr(fit, "selected_mediators")],
#'         WHO_indicator_codebook$indicator)]
#'
#' # A region-specific model (Europe).
#' des_eur <- WHO_mediation_design("imr", region = "EUR")
#'
#' @export
WHO_mediation_design <- function(outcome = c("imr", "u5mr", "leb",
                                             "lbw", "pou"),
                                 region = NULL, income = NULL,
                                 data = WHO_health_mediation) {
  outcome <- .WHO_match_arg(outcome, c("imr", "u5mr", "leb", "lbw", "pou"),
                            "outcome")
  .WHO_validate_data(data, outcome)
  df <- data
  if (!is.null(region)) {
    bad <- setdiff(region, levels(factor(df$region)))
    if (length(bad)) stop("Unknown region(s): ", paste(bad, collapse = ", "),
                          ".", call. = FALSE)
    df <- df[df$region %in% region, , drop = FALSE]
  }
  if (!is.null(income)) {
    # The groups `data` defines: a factor's levels (as before), or the
    # values of a character column (read.csv() and dplyr give one), which
    # has no levels and used to make every valid group "Unknown". Taken
    # from `data`, not from the region subset, so a region that lacks an
    # income group still reports "Too few complete observations (0)".
    known <- if (is.factor(data$income)) levels(data$income) else
      unique(as.character(data$income))
    bad <- setdiff(income, known)
    if (length(bad)) stop("Unknown income group(s): ", paste(bad, collapse = ", "),
                          ".", call. = FALSE)
    df <- df[as.character(df$income) %in% income, , drop = FALSE]
  }
  # Drop observations with no outcome or no exposure value.
  df <- df[!is.na(df[[outcome]]) & !is.na(df$gdp_growth), , drop = FALSE]
  if (nrow(df) < 10L)
    stop("Too few complete observations (", nrow(df),
         ") for this outcome/region/income combination.", call. = FALSE)
  .WHO_check_keys(df)

  indicators <- .WHO_indicator_names(data)
  X <- matrix(df$gdp_growth, ncol = 1L, dimnames = list(NULL, "gdp_growth"))
  Y <- df[[outcome]]
  M <- as.matrix(df[, indicators])

  # Confounders: region and income indicators (only if they still vary after
  # subsetting; unused factor levels dropped so no zero columns appear) plus
  # year as a numeric trend, the article's covariate set.
  terms <- character(0)
  reg <- droplevels(factor(as.character(df$region)))
  inc <- droplevels(factor(as.character(df$income)))
  if (nlevels(reg) > 1L) { df$.region <- reg; terms <- c(terms, ".region") }
  if (nlevels(inc) > 1L) { df$.income <- inc; terms <- c(terms, ".income") }
  terms <- c(terms, "year")
  Z <- stats::model.matrix(stats::reformulate(terms), data = df)[, -1,
                                                                 drop = FALSE]
  if (ncol(Z) == 0L) Z <- NULL

  list(X = X, Y = Y, M = M, Z = Z, outcome = outcome, n = nrow(df),
       n_countries = length(unique(df$code)),
       indicators = indicators, region = region, income = income)
}
