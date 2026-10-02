# Selected Mediators Under Each Multiplicity Method

Summarizes which mediators a
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
fit selected. When the fit was produced with
`report_all_methods = TRUE`, this reports the selected set under each of
the three multiplicity methods (Bonferroni for familywise error rate
control, Benjamini-Hochberg and Benjamini-Yekutieli for false discovery
rate control), so their conservativeness can be compared on the same
fit. Otherwise it reports the single method that was used.

## Usage

``` r
pe_selection(fit)
```

## Arguments

- fit:

  A result from
  [`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
  or one of its workers.

## Value

A tidy `data.frame` of class `poemed_tbl` (so its p-values print to
fixed decimals) with one row per multiplicity method and columns
`method`, `n_selected` (number of selected mediators), `j_pe` (the power
enhancement component \\J_m\\ under that method's screen), `stat_pe` and
`pval_pe` (the resulting power-enhanced statistic \\M\_{PE}\\ and its
p-value), and `selected_mediators` (their column names when `M` has
them, otherwise their column positions, comma-separated, or `"none"`).
Because each multiplicity method screens a different selected set into
\\J_m\\, the `pval_pe` column gives the PE / PE_BH / PE_BY global
p-values side by side, as the article's extended data-analysis tables
report them.

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
(its `report_all_methods` argument),
[`pe_mediators()`](https://yelleknek.github.io/POEMED/reference/pe_mediators.md).

Other mediation tests:
[`pe_mediate()`](https://yelleknek.github.io/POEMED/reference/pe_mediate.md),
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md),
[`pe_mediation_linear()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_linear.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_logistic.md),
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_poisson.md),
[`pe_mediators()`](https://yelleknek.github.io/POEMED/reference/pe_mediators.md),
[`summary.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/summary.poemed_tbl.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_mediation_data(n = 300, p = 100, outcome = "continuous",
                             pattern = "homogeneous", c1 = 1)
fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                    report_all_methods = TRUE)
#> Warning: HBIC selected lambda = 0.05, the smallest value of `lambda_grid`, for the full model. The criterion is minimized over the grid alone, so its minimum may lie beyond that end or between the end and its neighbor. Consider extending the grid past it and using a finer partition (more values), and compare the selected mediators and p-values across grids.
# On this homogeneous design (active mediators 1 to 5) the screens
# differ: Bonferroni, which controls the familywise error rate, selects
# three mediators; BH, the false discovery rate screen for independent or
# positively dependent mediators, selects five, one of them inactive; BY,
# valid under arbitrary dependence, selects the same three as Bonferroni.
pe_selection(fit)
#>  method     n_selected j_pe stat_pe pval_pe  selected_mediators
#>  Bonferroni 3          2332 2362    < 0.0001 3, 4, 5           
#>  BH         5          2954 2984    < 0.0001 2, 3, 4, 5, 14    
#>  BY         3          2332 2362    < 0.0001 3, 4, 5           
d$active_mediators
#> [1] 1 2 3 4 5
# HBIC chose 0.05, the smallest value of the default grid, so the fit
# warned. On a grid that extends below it and is spaced more finely,
# seq(0.01, 10, length.out = 200), it chooses about 0.11, and the
# screens then select 3 and 4 (Bonferroni) and 2, 3, and 4 (BH and BY):
# the selected set depends on the grid as well as on the screen.
```
