# WHO / World Bank Health-Expenditure Mediation Data (Benchmark Data Set)

A benchmark data set for high-dimensional mediation methods: the World
Health Organization (WHO) and World Bank panel analyzed in Yu and Kelley
(in press). It is a country-by-year panel used to ask how health-care
expenditure mediates the relationship between economic growth and
population health outcomes. The exposure is economic growth, the
candidate mediators are 57 health-expenditure indicators, and five
health outcomes are available. With 57 intercorrelated mediators this is
a genuinely high-dimensional mediation problem of the kind the
power-enhanced tests are built for, and the mediation effects are
heterogeneous (different kinds of spending push different outcomes in
different directions), which is exactly where the benchmark test loses
power and the PE test does not.

## Usage

``` r
WHO_health_mediation
```

## Format

`WHO_health_mediation` is a `data.frame` with 2002 rows (91 WHO member
states times the years 2000–2021) and 68 columns:

- code:

  ISO-3166 alpha-3 country code.

- country:

  Country name.

- region:

  WHO region: `AFR` (Africa, 32 members), `AMR` (Americas, 28), `EMR`
  (Eastern Mediterranean, 8), `EUR` (Europe, 9), `SEAR` (South-East
  Asia, 5), `WPR` (Western Pacific, 9).

- income:

  World Bank income group, a factor (not an ordered one) whose levels
  run from low to high: `Low` (16 members), `Lower-middle` (27),
  `Upper-middle` (26), `High` (22). The article abbreviates these Low,
  LM, UM, High.

- year:

  Calendar year, 2000–2021.

- gdp_growth:

  The exposure \\X\\: annual percentage growth of GDP per capita (World
  Bank Open Data), a common proxy for economic growth.

- imr:

  Infant mortality rate (deaths per 1000 live births).

- u5mr:

  Under-five mortality rate (deaths per 1000 live births).

- leb:

  Life expectancy at birth (years).

- lbw:

  Prevalence of low birthweight (percent of live births; covers 73
  members, 2000–2020, so other rows are `NA`).

- pou:

  Prevalence of undernourishment (percent of population; covers 79
  members, 2001–2021, so other rows are `NA`).

- che_gdp, ..., ext_usd2021_pc:

  The 57 candidate mediators \\M\\: health-expenditure indicators from
  the WHO Global Health Expenditure Database. Their codes and full names
  are in
  [WHO_indicator_codebook](https://yelleknek.github.io/POEMED/reference/WHO_indicator_codebook.md).
  They are complete (no missing values).

## Source

Exposure and health outcomes: World Bank Open Data, World Development
Indicators (<https://data.worldbank.org/>), licensed by the World Bank
under the Creative Commons Attribution 4.0 International license (CC BY
4.0, <https://datacatalog.worldbank.org/public-licenses>), which permits
redistribution and adaptation with attribution and a note of the changes
made. Health-expenditure mediators and the region and income labels:
World Health Organization, Global Health Expenditure Database
(<https://apps.who.int/nha/database>), made available under WHO's terms
and conditions for the use of its data
(<https://www.who.int/about/policies/publishing/data-policy/terms-and-conditions>):
non-commercial use for public health purposes, with attribution to WHO
and to the countries that provided the underlying data, and no
suggestion that WHO endorses any product or service. Attribution: World
Health Organization, Global Health Expenditure Database, accessed
November 2024; World Bank, World Development Indicators, accessed
November 2024. Changes made: the two sources were merged by country and
year, restricted to the 91 WHO member states with complete
health-expenditure data for 2000 to 2021, and the outcome series were
attached; no value was altered. The script that builds the data set,
`data-raw/WHO_health_mediation.R`, is in the package's GitHub repository
(<https://github.com/yelleKneK/POEMED>), not in the package itself. The
package's `LICENSE.note` file restates these terms.

## Details

It is shipped as a benchmark so that a method, old or new, can be
evaluated on a real, published high-dimensional mediation problem (see
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md),
whose output is the package's benchmark fit). The World Health
Organization (WHO) is the source of the health-expenditure mediators and
the region and income labels; the World Bank is the source of the
economic-growth exposure and the health outcomes.

Three of the five health outcomes (`imr`, `u5mr`, `leb`) cover all 91
members for 2000–2021; `lbw` and `pou` have narrower coverage, so they
carry `NA` where unavailable. The mediators are intercorrelated by
construction, since many are different normalizations of the same
underlying spending (per capita, percent of GDP, percent of current
health expenditure, constant versus current currency, and so on), which
is precisely the high-dimensional, correlated-mediator regime the
methods target.

Use
[`WHO_mediation_design()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_design.md)
to assemble the exposure, outcome, mediator, and confounder matrices for
a given outcome (and optionally a region or income subset) ready for
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md).
In the article's global model the confounders are the region, income
group, and year.

## References

Yu, X., & Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

## See also

[`WHO_mediation_design()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_design.md)
to build an analysis design,
[WHO_indicator_codebook](https://yelleknek.github.io/POEMED/reference/WHO_indicator_codebook.md)
for the mediator definitions,
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
to run the test.

Other WHO mediation data:
[`WHO_indicator_codebook`](https://yelleknek.github.io/POEMED/reference/WHO_indicator_codebook.md),
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md),
[`WHO_mediation_design()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_design.md)

## Examples

``` r
data(WHO_health_mediation)
table(WHO_health_mediation$region[WHO_health_mediation$year == 2010])
#> 
#>  AFR  AMR  EMR  EUR SEAR  WPR 
#>   32   28    8    9    5    9 

# A first test on the global infant-mortality design with the package
# defaults, which standardize X, M, and Z and center Y.
# WHO_mediation_analysis() instead applies the preprocessing and the
# tuning grid behind the article's tables, which leaves the
# confounders unscaled (see ?WHO_mediation_design for the difference).
# HBIC chooses 0.05, the smallest value of the default grid, so this fit
# warns; WHO_mediation_analysis() searches the article's grid instead.
des <- WHO_mediation_design("imr")
pe_mediation(des$X, des$Y, des$M, Z = des$Z, outcome = "continuous")
#> Warning: HBIC selected lambda = 0.05, the smallest value of `lambda_grid`, for the full model and the reduced model of the benchmark test. The criterion is minimized over the grid alone, so its minimum may lie beyond that end or between the end and its neighbor. Consider extending the grid past it and using a finer partition (more values), and compare the selected mediators and p-values across grids.
#>  term                  value    
#>  stat_hdmm             8.966e-06
#>  pval_hdmm             0.9976   
#>  stat_pe               238.6    
#>  j_pe                  238.6    
#>  pval_pe               < 0.0001 
#>  total_indirect_effect 0.0005369
#>  n_selected_mediators  2        
#>  df                    1        
#>  n_candidate_mediators 57       
#>  n_observations        2002     
#> 
#> Outcome model: continuous (linear)
#> Selected mediators (2): gge_gdp, ext_usd2021
#> Tuning parameter (HBIC): lambda = 0.05 from 100 values in [0.05, 10] (the grid's lower end)
```
