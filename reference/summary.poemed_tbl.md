# Summarize a POEMED Table

For a grouped table that compares the benchmark and power-enhanced tests
(any `poemed_tbl` with a `group` column, a `pval_hdmm` column, and a
power-enhanced p-value column (`pval_pe` or `pval_pe_bonferroni`), such
as the output of
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md)),
[`summary()`](https://rdrr.io/r/base/summary.html) digests it across
groups: how many groups each test detects mediation in, and in
particular the groups where the power-enhanced test detects mediation
that the benchmark misses, together with the most frequently flagged
mediators. For any other `poemed_tbl` it falls back to the ordinary
data-frame summary.

## Usage

``` r
# S3 method for class 'poemed_tbl'
summary(object, alpha_level = 0.05, ...)

# S3 method for class 'summary.poemed_tbl'
print(x, ...)
```

## Arguments

- object:

  A `poemed_tbl`.

- alpha_level:

  Significance level for counting a group as a detection. Default 0.05.

- ...:

  For a grouped comparison table, ignored. For any other `poemed_tbl`,
  passed to
  [`summary.data.frame()`](https://rdrr.io/r/base/summary.html) (for
  example `digits` or `quantile.type`).

- x:

  A `summary.poemed_tbl` object, for the print method.

## Value

For a grouped comparison table, an object of class `summary.poemed_tbl`
(a list, printed by its own method) with elements `outcome` (the outcome
label, or `NULL` when the table has none), `full_table` (logical; `TRUE`
for the extended layout of
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md)
with `full_table = TRUE`), `alpha_level`, `n_groups` (rows in the
table), `n_with_data` (groups with a benchmark p-value), `n_hdmm` and
`n_pe` (groups detected by each test), `pe_only` (a `data.frame` of the
groups the PE test detects but the benchmark does not, with their
`n_countries`, both p-values, and the `selected_mediators` codes), and
`top_mediators` (a frequency table of the mediators selected in the
groups the PE test detects). For any other `poemed_tbl`, the value of
[`summary.data.frame()`](https://rdrr.io/r/base/summary.html). The print
method returns its argument invisibly.

## Details

The printed digest follows the display convention of
[`print.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/poemed_tbl.md):
p-values appear to four decimal places, with a `< 0.0001` floor and
never in scientific notation, and the significance level appears as
supplied. The returned object keeps the p-values at full precision.

## See also

[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md),
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md).

Other mediation tests:
[`pe_mediate()`](https://yelleknek.github.io/POEMED/reference/pe_mediate.md),
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md),
[`pe_mediation_linear()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_linear.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_logistic.md),
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_poisson.md),
[`pe_mediators()`](https://yelleknek.github.io/POEMED/reference/pe_mediators.md),
[`pe_selection()`](https://yelleknek.github.io/POEMED/reference/pe_selection.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
# The global and income-group infant-mortality models, on a 25-point
# version of the default 100-point grid to keep the example fast (the
# global model's grid defaults to the same grid). The coarser grid
# selects the same mediators as the default, those of the article's
# Table 1, and moves some benchmark p-values slightly.
res <- WHO_mediation_analysis("imr", groupings = c("global", "income"),
                              lambda_grid = seq(0.1, 2, length.out = 25))
summary(res)
#> POEMED comparison: IMR
#>   groups: 5 (5 with data), alpha_level = 0.05
#>   detected by benchmark (HDMM): 1   by power-enhanced (PE): 3
#>   PE detects mediation in 2 groups the benchmark misses:
#>  group n_countries pval_hdmm  pval_pe selected_mediators
#>    ALL          91    0.0644 < 0.0001            gge_gdp
#>    Low          16    0.3089 < 0.0001       pvtd_usd2021
#>   most-flagged mediators:
#>     gge_gdp            1
#>     oops_che           1
#>     pvtd_usd2021       1
#>     shi_che            1
```
