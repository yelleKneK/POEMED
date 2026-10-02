# Broom Verbs for POEMED Results

`tidy()` and `glance()` methods so a POEMED result composes with the
broom ecosystem; the generics are re-exported, so `tidy(fit)` and
`glance(fit)` work after
[`library(POEMED)`](https://github.com/yelleKneK/POEMED) alone. For a
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
fit, `tidy()` returns the per-mediator table (see
[`pe_mediators()`](https://yelleknek.github.io/POEMED/reference/pe_mediators.md))
and `glance()` returns a one-row summary of the two tests and the total
indirect effect. For a
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md)
table both return the table itself (it is already one row per grid
point).

## Usage

``` r
# S3 method for class 'poemed_tbl'
tidy(x, ...)

# S3 method for class 'poemed_tbl'
glance(x, ...)

tidy(x, ...)

glance(x, ...)
```

## Arguments

- x:

  A `poemed_tbl` result.

- ...:

  Unused, for generic compatibility.

## Value

A `data.frame`. For a fit, `tidy()` has the columns of
[`pe_mediators()`](https://yelleknek.github.io/POEMED/reference/pe_mediators.md)
and `glance()` the columns `stat_hdmm`, `pval_hdmm`, `stat_pe`,
`pval_pe`, `total_indirect_effect` (one column per exposure when there
are several, suffixed with the exposure's column name or, for an unnamed
`X`, `_1`, `_2`, ...), `n_selected_mediators`, and `n_observations`.

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md),
[`pe_mediators()`](https://yelleknek.github.io/POEMED/reference/pe_mediators.md).

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_mediation_data(n = 120, p = 40, pattern = "contrasting",
                             outcome = "continuous", c1 = 1)
fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
#> Warning: HBIC selected lambda = 0.05, the smallest value of `lambda_grid`, for the full model and the reduced model of the benchmark test. The criterion is minimized over the grid alone, so its minimum may lie beyond that end or between the end and its neighbor. Consider extending the grid past it and using a finer partition (more values), and compare the selected mediators and p-values across grids.
# HBIC chose the smallest value of the default grid here, so the fit
# warned; see the lambda_grid argument of ?pe_mediation.
tidy(fit)
#>   mediator t_outcome t_exposure     screen_p selected
#> 1        1 15.551404  2.4315406 1.503477e-02    FALSE
#> 2        2 -6.550415  2.4345980 1.490834e-02    FALSE
#> 3        3  6.421123  5.1600037 2.469450e-07     TRUE
#> 4        4 -6.190102  5.2491552 1.527983e-07     TRUE
#> 5       35  1.546680 -0.5190176 6.037484e-01    FALSE
glance(fit)
#>   stat_hdmm pval_hdmm  stat_pe     pval_pe total_indirect_effect
#> 1  2.580174 0.1082097 417.6343 7.98664e-93             0.1270073
#>   n_selected_mediators n_observations
#> 1                    2            120
```
