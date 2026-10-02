#' WHO / World Bank Health-Expenditure Mediation Data (Benchmark Data Set)
#'
#' A benchmark data set for high-dimensional mediation methods: the World
#' Health Organization (WHO) and World Bank panel analyzed in Yu and
#' Kelley (in press). It is a country-by-year panel used to
#' ask how health-care expenditure mediates the relationship between
#' economic growth and population health outcomes. The exposure is
#' economic growth, the candidate mediators are 57 health-expenditure
#' indicators, and five health outcomes are available. With 57
#' intercorrelated mediators this is a genuinely high-dimensional
#' mediation problem of the kind the power-enhanced tests are built for,
#' and the mediation effects are heterogeneous (different kinds of
#' spending push different outcomes in different directions), which is
#' exactly where the benchmark test loses power and the PE test does not.
#'
#' It is shipped as a benchmark so that a method, old or new, can be
#' evaluated on a real, published high-dimensional mediation problem (see
#' [WHO_mediation_analysis()], whose output is the package's benchmark
#' fit). The World Health Organization (WHO) is the source
#' of the health-expenditure mediators and the region and income labels;
#' the World Bank is the source of the economic-growth exposure and the
#' health outcomes.
#'
#' @format `WHO_health_mediation` is a `data.frame` with 2002 rows (91
#'   WHO member states times the years 2000--2021) and 68 columns:
#' \describe{
#'   \item{code}{ISO-3166 alpha-3 country code.}
#'   \item{country}{Country name.}
#'   \item{region}{WHO region: `AFR` (Africa, 32 members), `AMR`
#'     (Americas, 28), `EMR` (Eastern Mediterranean, 8), `EUR` (Europe,
#'     9), `SEAR` (South-East Asia, 5), `WPR` (Western Pacific, 9).}
#'   \item{income}{World Bank income group, a factor (not an ordered one)
#'     whose levels run from low to high: `Low` (16 members),
#'     `Lower-middle` (27), `Upper-middle` (26), `High` (22). The article
#'     abbreviates these Low, LM, UM, High.}
#'   \item{year}{Calendar year, 2000--2021.}
#'   \item{gdp_growth}{The exposure \eqn{X}{X}: annual percentage growth of
#'     GDP per capita (World Bank Open Data), a common proxy for economic
#'     growth.}
#'   \item{imr}{Infant mortality rate (deaths per 1000 live births).}
#'   \item{u5mr}{Under-five mortality rate (deaths per 1000 live births).}
#'   \item{leb}{Life expectancy at birth (years).}
#'   \item{lbw}{Prevalence of low birthweight (percent of live births;
#'     covers 73 members, 2000--2020, so other rows are `NA`).}
#'   \item{pou}{Prevalence of undernourishment (percent of population;
#'     covers 79 members, 2001--2021, so other rows are `NA`).}
#'   \item{che_gdp, ..., ext_usd2021_pc}{The 57 candidate mediators
#'     \eqn{M}{M}: health-expenditure indicators from the WHO Global Health
#'     Expenditure Database. Their codes and full names are in
#'     [WHO_indicator_codebook]. They are complete (no missing values).}
#' }
#'
#' @details
#' Three of the five health outcomes (`imr`, `u5mr`, `leb`) cover all 91 members
#' for 2000--2021; `lbw` and `pou` have narrower coverage, so they carry
#' `NA` where unavailable. The mediators are intercorrelated by
#' construction, since many are different normalizations of the same
#' underlying spending (per capita, percent of GDP, percent of current
#' health expenditure, constant versus current currency, and so on), which
#' is precisely the high-dimensional, correlated-mediator regime the
#' methods target.
#'
#' Use [WHO_mediation_design()] to assemble the exposure, outcome,
#' mediator, and confounder matrices for a given outcome (and optionally a
#' region or income subset) ready for [pe_mediation()]. In the article's
#' global model the confounders are the region, income group, and year.
#'
#' @source Exposure and health outcomes: World Bank Open Data, World
#'   Development Indicators (\url{https://data.worldbank.org/}), licensed
#'   by the World Bank under the Creative Commons Attribution 4.0
#'   International license (CC BY 4.0,
#'   \url{https://datacatalog.worldbank.org/public-licenses}), which
#'   permits redistribution and adaptation with attribution and a note of
#'   the changes made. Health-expenditure mediators and the region and
#'   income labels: World Health Organization, Global Health Expenditure
#'   Database (\url{https://apps.who.int/nha/database}), made available
#'   under WHO's terms and conditions for the use of its data
#'   (\url{https://www.who.int/about/policies/publishing/data-policy/terms-and-conditions}):
#'   non-commercial use for public health purposes, with attribution to
#'   WHO and to the countries that provided the underlying data, and no
#'   suggestion that WHO endorses any product or service. Attribution:
#'   World Health Organization, Global Health Expenditure Database,
#'   accessed November 2024; World Bank, World Development Indicators,
#'   accessed November 2024. Changes made: the two sources
#'   were merged by country and year, restricted to the 91 WHO member
#'   states with complete health-expenditure data for 2000 to 2021, and
#'   the outcome series were attached; no value was altered. The script
#'   that builds the data set, `data-raw/WHO_health_mediation.R`, is in
#'   the package's GitHub repository
#'   (\url{https://github.com/yelleKneK/POEMED}), not in the package
#'   itself. The package's `LICENSE.note` file restates these terms.
#'
#' @references
#' Yu, X., & Kelley, K. (in press). Power Enhancement in
#' High-Dimensional Heterogeneous Mediation Analysis. \emph{Journal of the
#' American Statistical Association}.
#'
#' @seealso [WHO_mediation_design()] to build an analysis design,
#'   [WHO_indicator_codebook] for the mediator definitions,
#'   [pe_mediation()] to run the test.
#'
#' @family WHO mediation data
#'
#' @examples
#' data(WHO_health_mediation)
#' table(WHO_health_mediation$region[WHO_health_mediation$year == 2010])
#'
#' # A first test on the global infant-mortality design with the package
#' # defaults, which standardize X, M, and Z and center Y.
#' # WHO_mediation_analysis() instead applies the preprocessing and the
#' # tuning grid behind the article's tables, which leaves the
#' # confounders unscaled (see ?WHO_mediation_design for the difference).
#' # HBIC chooses 0.05, the smallest value of the default grid, so this fit
#' # warns; WHO_mediation_analysis() searches the article's grid instead.
#' des <- WHO_mediation_design("imr")
#' pe_mediation(des$X, des$Y, des$M, Z = des$Z, outcome = "continuous")
#'
#' @name WHO_health_mediation
#' @docType data
#' @keywords datasets
"WHO_health_mediation"

#' Codebook for the WHO Health-Expenditure Indicators
#'
#' The 57 health-expenditure indicators that serve as candidate mediators
#' in [WHO_health_mediation], with their short codes and descriptions.
#' The codes and their order are those of Table S.7 in the supplement to
#' Yu and Kelley (in press); the descriptions are adapted from that table
#' rather than copied, abbreviating some spelled-out denominators (for
#' example "as % of GDP" for "as % of Gross Domestic Product (GDP)" and
#' "as % of CHE" for "as % of Current Health Expenditure (CHE)").
#' Many indicators are alternative normalizations of the same underlying
#' spending (per capita, percent of GDP, percent of current health
#' expenditure, current versus constant 2021 currency, national currency
#' versus US dollars versus purchasing power parity), which is why they are
#' strongly intercorrelated.
#'
#' @format A `data.frame` with 57 rows and 2 columns:
#' \describe{
#'   \item{indicator}{The short indicator code, matching a column name in
#'     [WHO_health_mediation] (for example `gge_gdp`, `shi_che`).}
#'   \item{description}{The full indicator name. Common abbreviations:
#'     CHE (current health expenditure), GGHE-D (domestic general
#'     government health expenditure), PVT-D (domestic private health
#'     expenditure), EXT (external health expenditure), OOPS (out-of-pocket
#'     spending), SHI (social health insurance), GGE (general government
#'     expenditure), GDP (gross domestic product), NCU (national currency
#'     unit), PPP (purchasing power parity).}
#' }
#'
#' @source Adapted from Table S.7 ("A List of Health Spending Indicators
#'   and the Associated Indicator Code used in Real Data Analysis") of the
#'   supplement to Yu and Kelley (in press); underlying definitions from
#'   the WHO Global Health Expenditure Database.
#'
#' @references
#' Yu, X., & Kelley, K. (in press). Power Enhancement in
#' High-Dimensional Heterogeneous Mediation Analysis. \emph{Journal of the
#' American Statistical Association}.
#'
#' @seealso [WHO_health_mediation], [WHO_mediation_design()].
#'
#' @family WHO mediation data
#'
#' @examples
#' data(WHO_indicator_codebook)
#' head(WHO_indicator_codebook)
#' # Look up a mediator flagged in the article's analysis.
#' WHO_indicator_codebook[WHO_indicator_codebook$indicator == "gge_gdp", ]
#'
#' @name WHO_indicator_codebook
#' @docType data
#' @keywords datasets
"WHO_indicator_codebook"
