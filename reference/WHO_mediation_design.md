# Build a Mediation Design From the WHO Health-Expenditure Data

Assembles the exposure, outcome, mediator, and confounder matrices for
one of the analyses in Yu and Kelley (in press) from
[WHO_health_mediation](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md),
ready to pass to
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md).
The exposure is GDP per capita growth, the mediators are the 57
health-expenditure indicators, and the outcome is the chosen health
indicator. Optionally restrict to a WHO region and/or a World Bank
income group, as in the article's region-specific, income-specific, and
region-by-income models.

## Usage

``` r
WHO_mediation_design(
  outcome = c("imr", "u5mr", "leb", "lbw", "pou"),
  region = NULL,
  income = NULL,
  data = WHO_health_mediation
)
```

## Arguments

- outcome:

  Which health outcome to use as \\Y\\: `"imr"` (infant mortality rate),
  `"u5mr"` (under-five mortality rate), `"leb"` (life expectancy at
  birth), `"lbw"` (prevalence of low birthweight), or `"pou"`
  (prevalence of undernourishment).

- region:

  Optional character vector of WHO regions to keep (`"AFR"`, `"AMR"`,
  `"EMR"`, `"EUR"`, `"SEAR"`, `"WPR"`). Default `NULL` keeps all
  regions.

- income:

  Optional character vector of income groups to keep (`"Low"`,
  `"Lower-middle"`, `"Upper-middle"`, `"High"`). Default `NULL` keeps
  all groups.

- data:

  The source data frame, in the shape of
  [WHO_health_mediation](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md):
  the columns `code`, `region`, `income`, `year`, `gdp_growth`, and the
  chosen outcome, plus at least two indicator columns. Every column
  other than the key, exposure, and outcome columns of
  [WHO_health_mediation](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md)
  is treated as a candidate mediator, so the year, the exposure, the
  outcome, and every indicator column must be numeric; `region` and
  `income` may be character or factor. Every row the design uses needs
  its country `code` and `year`, and its `region` and `income` wherever
  those vary among the rows kept (that is, wherever they enter the
  confounders); a row missing one of them stops the design with a
  message naming the country. Rows outside the requested region or
  income group, or dropped for a missing outcome or exposure, are not
  checked. A `data` that breaks these rules stops with a message naming
  the column. Defaults to
  [WHO_health_mediation](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md).

## Value

A list with `X` (exposure, \\n \times 1\\), `Y` (outcome vector), `M`
(the \\n \times 57\\ mediator matrix), `Z` (the confounder matrix: the
region and income indicators that still vary after subsetting, and the
year), and the metadata `outcome`, `n` (country-year observations),
`n_countries` (unique countries), `indicators`, `region`, and `income`.

## Details

Rows with a missing outcome or missing exposure are dropped (the article
handles each outcome's coverage separately, which is why the
low-birthweight and undernourishment outcomes have fewer observations).
The confounder matrix \\Z\\ contains the WHO region and income group (as
indicator variables) and the year (numeric), the covariate set of the
article. The region and income indicators enter only while they still
vary after subsetting: a region-specific model drops the (now constant)
region, and unused factor levels are dropped so no all-zero indicator
columns remain. The year always enters, so a subset must span at least
two years; a single-year subset leaves a constant year column, which
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
rejects.

The indicator columns are passed through as they are. A missing value in
one of them makes
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
stop, since the method has no missing-data mechanism;
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md)
instead drops such an indicator within the affected group and names it
in a message.

To fit with the preprocessing of the analysis behind the article's
tables, standardize the exposure and mediators and center the outcome
while leaving the confounders (the region and income indicators and the
year) unscaled. The article's text says that covariates are
standardized; the unscaled confounders are a choice of that analysis.
The grid below, `seq(0.1, 2, length.out = 100)`, is the grid behind the
article's Tables 1, 2, and S.8 to S.12 and the default of
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md);
the call gives the first row of Table 1 (benchmark p-value 0.0644,
`gge_gdp` selected):


      des <- WHO_mediation_design("imr")
      pe_mediation(scale(des$X), des$Y - mean(des$Y), scale(des$M),
                   Z = des$Z, outcome = "continuous", scale = FALSE,
                   lambda_grid = seq(0.1, 2, length.out = 100))

HBIC chooses 0.1, the smallest value of that grid, for this fit, as it
does in about half of the article's cells, so the fit warns (class
`poemed_grid_boundary`);
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md)
keeps the article's grid and does not raise the warning (see its
Details).

The simpler call
`pe_mediation(des$X, des$Y, des$M, des$Z, outcome = "continuous")`
standardizes \\X\\, \\M\\, and \\Z\\ and centers \\Y\\, and so fits a
different model: scaling the confounder columns changes both the
benchmark Wald test and the penalized selection, so its p-values and its
selected set can differ from those of the call above.
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md)
applies the preprocessing and the grid of the article's analysis for
you.

## References

Yu, X., & Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md),
[WHO_health_mediation](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md),
[WHO_indicator_codebook](https://yelleknek.github.io/POEMED/reference/WHO_indicator_codebook.md).

Other WHO mediation data:
[`WHO_health_mediation`](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md),
[`WHO_indicator_codebook`](https://yelleknek.github.io/POEMED/reference/WHO_indicator_codebook.md),
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
# Global infant-mortality model: GDP growth -> health spending -> IMR,
# with the preprocessing and the grid of the article's Table 1; the fit
# selects gge_gdp. HBIC chooses 0.1, the smallest value of that grid, so
# the fit warns (class poemed_grid_boundary); WHO_mediation_analysis()
# keeps the article's grid and does not raise the warning.
des <- WHO_mediation_design("imr")
c(n = des$n, mediators = ncol(des$M))
#>         n mediators 
#>      2002        57 
fit <- pe_mediation(scale(des$X), des$Y - mean(des$Y), scale(des$M),
  Z = des$Z, outcome = "continuous", scale = FALSE,
  lambda_grid = seq(0.1, 2, length.out = 100))
#> Warning: HBIC selected lambda = 0.1, the smallest value of `lambda_grid`, for the full model and the reduced model of the benchmark test. The criterion is minimized over the grid alone, so its minimum may lie beyond that end or between the end and its neighbor. Consider extending the grid past it and using a finer partition (more values), and compare the selected mediators and p-values across grids.
fit
#>  term                  value   
#>  stat_hdmm             3.421   
#>  pval_hdmm             0.0644  
#>  stat_pe               145     
#>  j_pe                  141.6   
#>  pval_pe               < 0.0001
#>  total_indirect_effect 0.372   
#>  n_selected_mediators  1       
#>  df                    1       
#>  n_candidate_mediators 57      
#>  n_observations        2002    
#> 
#> Outcome model: continuous (linear)
#> Selected mediators (1): gge_gdp
#> Tuning parameter (HBIC): lambda = 0.1 from 100 values in [0.1, 2] (the grid's lower end)
# The selected mediator's full name:
WHO_indicator_codebook$description[
  match(des$indicators[attr(fit, "selected_mediators")],
        WHO_indicator_codebook$indicator)]
#> [1] "General Government Expenditure (GGE) as % of GDP"

# A region-specific model (Europe).
des_eur <- WHO_mediation_design("imr", region = "EUR")
```
