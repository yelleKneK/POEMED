# Run the Article's Empirical Mediation Analysis

Runs the article's real-data analysis for one health outcome across the
groupings it reports: a global model using all WHO members, then
separate models within each WHO region, within each World Bank income
group, and (optionally) within each region-by-income cell. For every
grouping it reports the benchmark Wald test and the power-enhanced test
on the total indirect effect of health expenditure, together with the
individual mediators the PE screen selects. The result has the layout of
the article's data-analysis tables (Table 1 for infant mortality and
Table 2 for undernourishment; Tables S.8 to S.12 of its supplement give
the extended layout and the other three outcomes), computed from the
shipped
[WHO_health_mediation](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md)
data.

## Usage

``` r
WHO_mediation_analysis(
  outcome = c("imr", "u5mr", "leb", "lbw", "pou"),
  groupings = c("global", "region", "income"),
  error_level = 0.05,
  lambda_grid = seq(0.1, 2, length.out = 100),
  lambda_grid_global = lambda_grid,
  data = WHO_health_mediation,
  cores = 1L,
  verbose = FALSE,
  full_table = FALSE
)
```

## Arguments

- outcome:

  Which health outcome to analyze: `"imr"`, `"u5mr"`, `"leb"`, `"lbw"`,
  or `"pou"`. See
  [WHO_health_mediation](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md).

- groupings:

  Character vector of groupings to run, any of `"global"`, `"region"`,
  `"income"`, and `"region_income"`. Default
  `c("global", "region", "income")` (the region-by-income cells are the
  slowest and are opt-in).

- error_level:

  Target familywise error rate for the mediator-screening step. Default
  0.05, as in Yu and Kelley (in press).

- lambda_grid:

  Tuning-parameter grid for the subgroup penalized fits (regional,
  income, and region-by-income models). Default
  `seq(0.1, 2, length.out = 100)`, the grid behind the article's
  real-data tables (Tables 1, 2, and S.8 to S.12 of Yu and Kelley, in
  press); it differs from the default of
  [`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
  because the function reproduces the article's analysis as run.

- lambda_grid_global:

  Tuning-parameter grid for the global (all-country) model only.
  Defaults to `lambda_grid`: the article's tables use the same grid for
  the global model and the subgroup models.

- data:

  The source data, in the shape of
  [WHO_health_mediation](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md);
  see
  [`WHO_mediation_design()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_design.md)
  for the columns it must carry. Defaults to
  [WHO_health_mediation](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md).

- cores:

  Number of CPU cores; values above 1 fit the subgroups in parallel via
  the base parallel package (Unix only). Default 1.

- verbose:

  Logical; if `TRUE`, print each group as it is fit. Default `FALSE`.

- full_table:

  Logical; if `TRUE`, return the article's *extended* table layout, with
  the power-enhanced p-value and the selected mediator codes reported
  separately under each multiplicity method (Bonferroni for FWER, BH and
  BY for FDR) instead of a single primary method. Default `FALSE`.

## Value

A `data.frame` with one row per fitted subgroup. With
`full_table = FALSE` (the default) the columns are `grouping` (global /
region / income / region_income), `group` (the specific group, for
example `AFR` or `Low`), `n_countries` (number of countries; `NA` for a
group skipped at the design stage: one with fewer than 10 complete
observations, which includes every group the article's tables list with
0 countries, or one whose rows lack a key the group needs), `pval_hdmm`
and `pval_pe` (the benchmark and power-enhanced p-values for the total
indirect effect), `n_selected` (number of mediators the screen
selected), and `selected_mediators` (their codes, or `"none"`). With
`full_table = TRUE` the single
`pval_pe`/`n_selected`/`selected_mediators` columns are replaced by
`pval_pe_bonferroni`, `pval_pe_bh`, `pval_pe_by` and the matching
`selected_bonferroni`, `selected_bh`, `selected_by` code columns. The
result is a `poemed_tbl` with the attributes `"outcome"`,
`"full_table"`, and `"notes"` (see Details); call
[summary()](https://yelleknek.github.io/POEMED/reference/summary.poemed_tbl.md)
on it for a cross-group digest of where the power-enhanced test detects
mediation the benchmark misses.

## Details

The output is the package's **benchmark fit**: the test suite pins the
selected sets and the benchmark p-values of the article's tables under
the function's default grid, so a change to the estimation code that
moved them would be caught, and the fit serves as a fixed reference
point for the method.

Each subgroup is fit with the preprocessing of the analysis behind the
article's tables: the exposure and the 57 mediators are standardized,
the outcome is centered, and the confounders (the region and income
indicators and the year) are left unscaled. The article's text says that
covariates are standardized; the unscaled confounders are a choice of
that analysis. Within a group, an indicator that is constant there
(common in the small region-by-income cells) is dropped before fitting,
since it carries no information. An indicator with a missing or
non-finite value among the group's rows is dropped as well, because the
method has no missing-data mechanism; those are named in a message,
since they change the set of candidate mediators.

Four kinds of event are reported once, after all groups are fit, and
recorded in the `"notes"` attribute (a `data.frame` with `group`,
`stage`, and `message`, or `NULL` when there are none). The `stage`
column says which:

- `"design"`: a group with too few complete observations (some
  region-by-income cells are empty or nearly so), or whose rows lack a
  country code, a year, or a region or income group that varies within
  the group, is skipped and reported as an `NA` row, named in a message.

- `"missing"`: indicators with missing or non-finite values were dropped
  in the group before fitting, named in a message.

- `"empty"`: the penalized fit that HBIC selected contains no mediator.
  The row keeps the values of POEMED's reporting convention for an empty
  selection (p-values of 1, no selected mediator), which are not
  evidence for the null, and the group is named, with the lambda HBIC
  selected, in a warning of class `poemed_empty_fit` at any value of
  `cores`.

- `"fit"`: the fit failed; the group is reported as an `NA` row and
  named in a warning that quotes the error.

Because every group fits a penalized model over a grid of tuning
parameters, a full run with all groupings fits dozens of models and
takes about ten seconds on a current laptop; start with the default
groupings, or a single outcome, before scaling up. In 67 of the 128
cells that the five outcomes fit, the HBIC choice for the full penalized
model or for the reduced model of the benchmark Wald test (54 cells
each) is 0.1, the grid's smallest value; that is a property of the
article's analysis on its fixed grid, so the function does not raise the
single-fit `poemed_grid_boundary` warning for it. To see a cell's tuning
record, fit the cell by hand with
[`pe_mediation_linear()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_linear.md),
as the reproduction vignette shows.

The tuning-parameter grid is `seq(0.1, 2, length.out = 100)` for the
global model and for every subgroup model, the grid behind the article's
Tables 1, 2, and S.8 to S.12 (Section 5 of Yu and Kelley, in press); the
article's text itself states only that the tuning parameter is chosen by
the HBIC criterion. The article's real-data analysis was run with this
package
([`pe_mediation_linear()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_linear.md)
at that grid, with the preprocessing described above);
`WHO_mediation_analysis()` automates those calls. Under this default the
function reproduces those tables. The selected set agrees with the
printed one in every fitted cell under Bonferroni, BY, and BH, with one
exception (the BH set for under-five mortality in the High income group,
where the article's set also holds `gge_gdp`). The benchmark p-value
agrees to four decimals in all but a handful of cells: infant mortality
in SEAR (0.6483 here, 0.6482 printed) and in AMR-High (0.0726 here,
0.0725 printed) differ by rounding in the fourth decimal, and a few
p-values at or below 0.0001, which the article prints to one significant
figure, come out somewhat smaller here. Exact p-values depend on the
grid and on the unscaled confounders. Six of the region-by-income cells
that the article's tables fill are not fit at all and are reported as
`NA` (the cells the tables list with 0 countries are `NA` too, skipped
before fitting):

- four cells in which two indicators are exact linear partners and the
  fit keeps both, so a matrix it inverts is singular: infant mortality
  and under-five mortality in EMR-Lower-middle (`dom_che` and `ext_che`
  sum to 100 there), and life expectancy and under-five mortality in
  AFR-Upper-middle (`ext_che` with `dom_che`, and `shi_che` with
  `chi_che`). The article's values are reproduced by dropping one
  partner; see the recipe below.

- two cells of the undernourishment outcome, EUR-High and WPR-High, in
  which the outcome is 2.5 in every country-year (the lower reporting
  bound of the source), so there is no variation to mediate. The
  function reports the fit as failed with the message that `Y` is
  constant; the article prints a p-value of 1 there.

Many of the 57 indicators are exact or near-exact linear combinations of
one another, and the penalized selection chooses among such partners. As
a result the benchmark p-values, and occasionally the selected
mediators, can change when the indicator columns of `data` are
reordered, even between two indicators that are exactly collinear within
a group. The results reported are those for the column order of the
shipped data.

A group whose fit fails with a singular-matrix error usually holds two
indicators that are exact linear partners within that group (for
example, `dom_che` and `ext_che` sum to 100 in the Eastern Mediterranean
lower-middle-income group;
[`stats::cor()`](https://rdrr.io/r/stats/cor.html) of 1 or -1 among the
columns of the group's `M` finds such pairs). Such a group can be fit by
hand after dropping one member of the pair. Which member is dropped
matters, since the two fits select among different columns (dropping
`ext_che` instead in the infant-mortality cell below gives 0.0786). The
article's analysis removed `dom_che` and `vfa_che` in the two
EMR-Lower-middle cells (`cfa_che` and `vfa_che` are near-exact partners
there as well), giving benchmark p-values of 0.0837 for infant mortality
and 0.0914 for under-five mortality with no mediator selected, and
`shi_che` and `ext_che` in the two AFR-Upper-middle cells, giving 0.1311
for life expectancy with no mediator selected and 0.0389 for under-five
mortality with `pvtd_ncu2021_pc` selected. With the default grid, the
infant-mortality cell is:


      cell <- WHO_mediation_design("imr", region = "EMR",
                                   income = "Lower-middle")
      keep <- apply(cell$M, 2, sd) > 0 &
        !(cell$indicators 
      grid <- seq(0.1, 2, length.out = 100)
      pe_mediation(scale(cell$X), cell$Y - mean(cell$Y),
                   scale(cell$M[, keep]), Z = cell$Z,
                   outcome = "continuous", scale = FALSE,
                   lambda_grid = grid)

which prints the benchmark p-value 0.0837 of the article's Table 1.

## How to Cite

If you use POEMED in published work, please cite Yu and Kelley (in
press), the article that introduces its methods and the package.
`citation("POEMED")` gives the full reference and a BibTeX entry.

## References

Yu, X., & Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

## See also

[`WHO_mediation_design()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_design.md)
for a single design,
[WHO_health_mediation](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md)
for the data,
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
for the test.

Other WHO mediation data:
[`WHO_health_mediation`](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md),
[`WHO_indicator_codebook`](https://yelleknek.github.io/POEMED/reference/WHO_indicator_codebook.md),
[`WHO_mediation_design()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_design.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
# The global infant-mortality model, the first row of the article's
# Table 1 (under a second): benchmark p-value 0.0644, power-enhanced
# p-value below 0.0001, and gge_gdp selected. Add "region" and
# "income" to groupings for the rest of the table.
WHO_mediation_analysis("imr", groupings = "global")
#>  grouping group n_countries pval_hdmm pval_pe  n_selected selected_mediators
#>  global   ALL   91          0.0644    < 0.0001 1          gge_gdp           
#> 
#> Outcome: imr
```
