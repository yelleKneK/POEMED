# Per-Mediator Detail Behind a Power-Enhanced Mediation Test

Returns the per-mediator table that explains a
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
result. The table has one row for each candidate mediator, that is, each
mediator the penalized fit retained with a nonzero coefficient. Each row
gives the mediator's two path statistics, its screening p-value, and
whether the screen selected it. This is the "why" behind the selected
set, so a user can see which candidates were close to the screening
threshold and which were screened out.

## Usage

``` r
pe_mediators(fit)
```

## Arguments

- fit:

  A result from
  [`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
  or one of its workers.

## Value

A tidy `data.frame` (class `poemed_tbl`) with columns:

- mediator:

  Column position of the mediator in the original `M`.

- name:

  The mediator's column name in `M`; present only when `M` has column
  names.

- t_outcome:

  The standardized mediator-on-outcome statistic (the \\M \to Y\\ path,
  \\\hat\alpha\_{m,j} / \hat\sigma\_{m,j}\\).

- t_exposure:

  The standardized exposure-on-mediator statistic (the \\X \to M\\
  path). With one exposure it is the signed *t* statistic, so a negative
  value marks a negative path. With several exposures it is the largest
  absolute *t* statistic across the exposures, so it is never negative
  and does not show the path's direction; regress the mediator on the
  exposures and any confounders (for example with
  [`lm()`](https://rdrr.io/r/stats/lm.html)) to see the signs.

- screen_p:

  The screening p-value, the larger of the two path p-values (the
  smaller of these across exposures when there are several). A mediator
  is selected when this clears the multiplicity threshold.

- selected:

  Logical; `TRUE` when the screen selected the mediator under the fit's
  `method`, so the `TRUE` rows are the fit's selected set. Every row is
  a candidate the penalized fit retained, so `FALSE` marks a candidate
  the screen did not select.

The table is empty when the penalized fit retained no candidate
mediators.

## Details

Two steps narrow the mediators, and the table keeps them apart. The
penalized fit keeps the candidates, which are the rows of the table. The
multiplicity screen then selects mediators among those candidates, which
are the rows whose `selected` column is `TRUE`. A row with
`selected = FALSE` is therefore a candidate that the screen did not
select, not a mediator the penalized fit dropped.

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md),
[`pe_selection()`](https://yelleknek.github.io/POEMED/reference/pe_selection.md).

Other mediation tests:
[`pe_mediate()`](https://yelleknek.github.io/POEMED/reference/pe_mediate.md),
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md),
[`pe_mediation_linear()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_linear.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_logistic.md),
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_poisson.md),
[`pe_selection()`](https://yelleknek.github.io/POEMED/reference/pe_selection.md),
[`summary.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/summary.poemed_tbl.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                             pattern = "contrasting", c1 = 1)
fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
#> Warning: HBIC selected lambda = 0.05, the smallest value of `lambda_grid`, for the full model and the reduced model of the benchmark test. The criterion is minimized over the grid alone, so its minimum may lie beyond that end or between the end and its neighbor. Consider extending the grid past it and using a finer partition (more values), and compare the selected mediators and p-values across grids.
# HBIC chose the smallest value of the default grid here, so the fit
# warned; see the lambda_grid argument of ?pe_mediation.
pe_mediators(fit)
#>  mediator t_outcome t_exposure screen_p selected
#>  1        21.75     0.3641     0.7158   FALSE   
#>  2        -10.21    3.131      0.0017    TRUE   
#>  3        8.602     4.349      < 0.0001  TRUE   
#>  4        -6.82     6.185      < 0.0001  TRUE   
#>  33       -2.378    -13.08     0.0174   FALSE   
```
