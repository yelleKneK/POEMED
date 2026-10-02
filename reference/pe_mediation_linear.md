# Power-Enhanced Mediation Test for a Continuous Outcome

Tests the global null hypothesis of no active mediator among a
high-dimensional set of candidate mediators, for a continuous
(linear-model) outcome, and reports both the benchmark Wald test on the
total indirect effect and the power-enhanced (PE) test of Yu and Kelley
(in press), which remains powerful when individual mediation effects are
heterogeneous or contrasting. This is the continuous-outcome worker
behind
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md);
call it directly when the outcome is continuous.

## Usage

``` r
pe_mediation_linear(
  X,
  Y,
  M,
  Z = NULL,
  method = c("Bonferroni", "BH", "BY"),
  scale = TRUE,
  error_level = 0.05,
  report_all_methods = FALSE,
  drop_constant = FALSE,
  lambda_grid = seq(0.05, 10, length.out = 100)
)
```

## Arguments

- X:

  Numeric matrix of exposures with \\n\\ rows and \\q\\ columns (\\q\\
  is usually 1). A numeric vector is treated as a single exposure.

- Y:

  Numeric outcome vector of length \\n\\. Continuous for
  `outcome = "continuous"`, coded 0/1 for `"binary"`, and nonnegative
  integer counts for `"count"`.

- M:

  Numeric matrix of candidate mediators with \\n\\ rows and \\p\\
  columns; \\p\\ may exceed \\n\\. Column names, when present, label the
  selected mediators in the print footer and the mediator table; without
  them mediators are reported by column position.

- Z:

  Optional numeric matrix of confounders with \\n\\ rows. Default `NULL`
  (no confounders).

- method:

  Multiplicity adjustment used to screen individual mediators in the PE
  component: `"Bonferroni"` (familywise error rate control, the
  default), `"BH"`, or `"BY"` (false discovery rate control by the
  procedures of Benjamini and Hochberg, 1995, and of Benjamini and
  Yekutieli, 2001; `"BH"` assumes independence or positive dependence
  among mediators, `"BY"` is valid under arbitrary dependence). See
  Details for the screening threshold each one uses.

- scale:

  Logical; if `TRUE` (default) the exposures, mediators, and confounders
  are each standardized to mean 0 and standard deviation 1 before
  fitting, as in the article, and a continuous outcome is centered but
  not rescaled. The total indirect effect is then the change in the
  outcome per standard deviation of each exposure: in the units of `Y`
  for a continuous outcome, in log odds for a binary outcome, and in log
  mean count for a count outcome. With `FALSE` the inputs are used
  exactly as supplied, and the total indirect effect is per unit of each
  exposure. The continuous-outcome fit has no intercept, and the default
  tuning grid was chosen for standardized columns, so `FALSE` is meant
  for inputs that are already standardized (and, for a continuous
  outcome, a `Y` that is already centered); on raw inputs it fits a
  different model and can even change the sign of the estimate. Because
  a continuous `Y` is not rescaled, the tuning grid is not expressed in
  the outcome's units, so the selected mediators, both p-values, and the
  verdict depend on the units of `Y`: the same data with `Y` multiplied
  by 10 can select a different set of mediators.

- error_level:

  Target error rate for the mediator-screening step, in \\(0, 1)\\.
  Interpreted as the familywise error rate when `method = "Bonferroni"`
  and as the false discovery rate when `method = "BH"` or `"BY"`.
  Default 0.05. This is distinct from the significance level used to
  test the global null (which is the user's choice when reading
  `pval_pe`). The PE test is asymptotically valid for any fixed value in
  \\(0, 1)\\, but in finite samples a larger value lets the screen pass
  more spurious mediators under the null and raises the PE test's size
  further above the benchmark's (see Details).

- report_all_methods:

  Logical; if `TRUE`, the selected-mediator set is computed for all
  three multiplicity methods (Bonferroni, BH, BY) off the same fit and
  recorded for comparison, retrievable with
  [`pe_selection()`](https://yelleknek.github.io/POEMED/reference/pe_selection.md).
  The `method` argument still drives the primary result. Default
  `FALSE`.

- drop_constant:

  Logical; how to handle constant (zero-variance) mediator columns,
  which cannot be standardized or carry signal. If `FALSE` (the default)
  they are an error; if `TRUE` they are dropped with a warning and the
  remaining mediators are renumbered back to their original column
  positions in the reported selected set. Useful for small subgroups
  where some mediators happen to be constant.

- lambda_grid:

  Numeric vector of candidate tuning parameters for the SCAD penalty
  (Fan and Li, 2001) in the penalized mediator fit, which for a binary
  outcome is adaptively rescaled (see
  [`pe_mediation_logistic()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_logistic.md)).
  The fit is repeated at each value and the one minimizing the
  high-dimensional BIC (HBIC) of Wang, Kim, and Li (2013) is kept:
  whatever the grid, the value selected is the best of the values
  offered, since HBIC is minimized over the grid and not over every
  positive lambda. The default, `seq(0.05, 10, length.out = 100)` for
  every outcome family, is written out in the signature so that it can
  be edited in place. Every value is fit (the search is not adaptive),
  so the cost of a fit grows with the length of the grid. The grid can
  be tuned to the situation. Fewer values, or a range narrowed to where
  earlier fits put the HBIC minimum, make each fit faster, which matters
  when many fits are run (a simulation study, for example); a wider or
  finer grid is worth trying when the selected value is the smallest or
  largest value of the grid, in which case a non-empty fit warns (class
  `poemed_grid_boundary`) and the print footer, which reports the full
  model's choice, says so. For a continuous outcome the benchmark Wald
  test also needs a penalized fit of the reduced model, the outcome on
  the mediators and confounders without the exposure; that fit searches
  the same grid, its choice is `lambda_selected_reduced` in the
  `"tuning"` attribute, and the warning covers it as well, naming
  whether the full model, the reduced model, or both sat at the end. For
  a continuous outcome the grid is fixed while `Y` is centered but not
  rescaled, so the value HBIC chooses, and with it the selected
  mediators, depends on the units of `Y` (see `scale`). The grid
  searched and the value chosen are reported in the `"tuning"` attribute
  and the print footer.

## Value

A tidy `data.frame` of class `poemed_tbl` with rows `stat_hdmm` and
`pval_hdmm` (the benchmark Wald test), `stat_pe`, `j_pe`, and `pval_pe`
(the power-enhanced test and its PE component), `total_indirect_effect`
(the estimated total indirect effect; one row per exposure when `q > 1`,
suffixed with the exposure's column name when `X` has column names and
with `_1`, `_2`, ... otherwise), `n_selected_mediators` (the number of
mediators the screen selected), `df` (the chi-square degrees of freedom,
equal to `q`), `n_candidate_mediators`, and `n_observations`. With
`scale = TRUE` the total indirect effect is in the units of `Y` per
standard deviation of the exposure; with `scale = FALSE` it is in the
units of `Y` per unit of the exposure as supplied. No confidence
interval is reported: the article reports none, and the Wald interval
built on the penalized estimate under-covers in finite samples.

Attributes: `"selected_mediators"` (the column positions in `M` that the
screen selected; the print footer shows their column names when `M` has
them), `"mediator_names"` and `"exposure_names"` (the column names of
`M` and `X`, or `NULL` when unnamed), `"mediator_table"` (the
per-mediator screening statistics for the penalized-selected mediators,
see
[`pe_mediators()`](https://yelleknek.github.io/POEMED/reference/pe_mediators.md);
its `t_exposure` is the signed exposure-on-mediator statistic when
`q = 1` and the largest absolute one across exposures when `q > 1`, and
its `screen_p` is the smallest screening p-value across exposures),
`"method"`, `"error_level"`, `"outcome"`, `"p_terms"` (the names of the
rows that hold p-values, which the print method formats as p-values),
`"tuning"` (a list with the grid searched, the HBIC-selected
`lambda_selected`, and `at_lower_end` and `at_upper_end`, `TRUE` when
the selection sat at the grid's smallest or largest value; for a
continuous outcome also `lambda_selected_reduced`, the choice for the
reduced model of the benchmark Wald test, searched on the same grid),
`"empty_fit"`, and `"pe_by_method"` (a data frame with one row per
multiplicity method screened, giving that method's `stat_pe`, `j_pe`,
`pval_pe`, and number of selected mediators; it has one row unless
`report_all_methods = TRUE`). With `report_all_methods = TRUE` there is
also `"selection_by_method"`, the selected set under each method (see
[`pe_selection()`](https://yelleknek.github.io/POEMED/reference/pe_selection.md)).

When the penalized fit that HBIC selects contains no mediator there is
nothing to test: the function warns (a condition of class
`poemed_empty_fit`), sets `"empty_fit"` to `TRUE`, and returns the same
rows with both statistics 0, both p-values 1 (POEMED's reporting
convention, under which an empty selection is a non-rejection; the
article states only that \\J_m = 0\\ when the selection is empty), and a
total indirect effect of 0. A p-value of 1 from such a fit is not
evidence for the null. An empty fit is expected when no mediator is
active. A non-empty fit warns (class `poemed_grid_boundary`) when its
HBIC-selected lambda is the smallest or largest value of the grid, for
the full model or, for a continuous outcome, for the reduced model of
the benchmark Wald test, and the warning names the fit that sat there:
the criterion is minimized over the grid alone, so extend the grid past
that end and space its values more finely.

## Details

The model is the linear mediation pair \$\$Y = \alpha_m' M + \alpha_x'
X + \alpha_z' Z + \varepsilon_y, \qquad M = \Gamma_x' X + \Gamma_z' Z +
\varepsilon_m,\$\$ with total indirect effect \\\beta = \Gamma_x
\alpha_m\\. Following Guo et al. (2023), the mediator coefficients
\\\alpha_m\\ are estimated by partial penalized least squares with a
SCAD penalty (Fan and Li, 2001), leaving the exposures and confounders
unpenalized. The penalized fit is computed by the one-step local linear
approximation of Zou and Li (2008) started from a lasso fit, as in Fan,
Xue, and Zou (2014): a lasso fit, then a weighted lasso fit whose
weights are the SCAD derivatives at the lasso coefficients. The tuning
parameter is chosen by the high-dimensional BIC of Wang, Kim, and Li
(2013). The fit has no intercept, because with `scale = TRUE` the
outcome is centered and the other inputs standardized first (see `scale`
in
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)).

The benchmark statistic is the Wald statistic \\S_n = n \hat\beta'
\hat\Sigma\_\beta^{-1} \hat\beta\\ of Guo et al. (2022), which extends
the partially penalized Wald test of Guo et al. (2023) to observed
confounders. The power-enhanced statistic of Yu and Kelley (in press),
built on the power enhancement principle of Fan, Liao, and Yao (2015),
adds the component \\J_m\\ that accumulates the marginal signal from
each selected mediator (see
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
for the rationale, the formula, the screening thresholds, and the test's
finite-sample size). Both are referred to a chi-square distribution with
\\q\\ (the number of exposures) degrees of freedom, their common limit
under the global null.

## How to Cite

If you use POEMED in published work, please cite Yu and Kelley (in
press), the article that introduces its methods and the package.
`citation("POEMED")` gives the full reference and a BibTeX entry.

## References

Yu, X., & Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

Guo, X., Li, R., Liu, J., & Zeng, M. (2022). High-dimensional mediation
analysis for selecting DNA methylation loci mediating childhood trauma
and cortisol stress reactivity. *Journal of the American Statistical
Association, 117*(539), 1110–1121.
[doi:10.1080/01621459.2022.2053136](https://doi.org/10.1080/01621459.2022.2053136)

Guo, X., Li, R., Liu, J., & Zeng, M. (2023). Statistical inference for
linear mediation models with high-dimensional mediators and application
to studying stock reaction to COVID-19 pandemic. *Journal of
Econometrics, 235*(1), 166–179.
[doi:10.1016/j.jeconom.2022.03.001](https://doi.org/10.1016/j.jeconom.2022.03.001)

Fan, J., Liao, Y., & Yao, J. (2015). Power enhancement in
high-dimensional cross-sectional tests. *Econometrica, 83*(4),
1497–1541. [doi:10.3982/ECTA12749](https://doi.org/10.3982/ECTA12749)

Fan, J., & Li, R. (2001). Variable selection via nonconcave penalized
likelihood and its oracle properties. *Journal of the American
Statistical Association, 96*(456), 1348–1360.
[doi:10.1198/016214501753382273](https://doi.org/10.1198/016214501753382273)

Zou, H., & Li, R. (2008). One-step sparse estimates in nonconcave
penalized likelihood models. *The Annals of Statistics, 36*(4),
1509–1533.
[doi:10.1214/009053607000000802](https://doi.org/10.1214/009053607000000802)

Fan, J., Xue, L., & Zou, H. (2014). Strong oracle optimality of folded
concave penalized estimation. *The Annals of Statistics, 42*(3),
819–849. [doi:10.1214/13-AOS1198](https://doi.org/10.1214/13-AOS1198)

Wang, L., Kim, Y., & Li, R. (2013). Calibrating nonconvex penalized
regression in ultra-high dimension. *The Annals of Statistics, 41*(5),
2505–2536. [doi:10.1214/13-AOS1159](https://doi.org/10.1214/13-AOS1159)

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
rate: A practical and powerful approach to multiple testing. *Journal of
the Royal Statistical Society, Series B, 57*(1), 289–300.
[doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

Benjamini, Y., & Yekutieli, D. (2001). The control of the false
discovery rate in multiple testing under dependency. *The Annals of
Statistics, 29*(4), 1165–1188.
[doi:10.1214/aos/1013699998](https://doi.org/10.1214/aos/1013699998)

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
for the front end,
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_logistic.md)
and
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_poisson.md)
for binary and count outcomes, and
[`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md)
to generate data for trying the method.

Other mediation tests:
[`pe_mediate()`](https://yelleknek.github.io/POEMED/reference/pe_mediate.md),
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_logistic.md),
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_poisson.md),
[`pe_mediators()`](https://yelleknek.github.io/POEMED/reference/pe_mediators.md),
[`pe_selection()`](https://yelleknek.github.io/POEMED/reference/pe_selection.md),
[`summary.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/summary.poemed_tbl.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
# A contrasting setting: four active mediators whose indirect effects
# cancel, so the total indirect effect is zero. The Wald test does not
# reject here; the PE test rejects and selects three of the four
# active mediators.
d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
                             outcome = "continuous", c1 = 1)
pe_mediation_linear(d$X, d$Y, d$M)
#> Warning: HBIC selected lambda = 0.05, the smallest value of `lambda_grid`, for the full model and the reduced model of the benchmark test. The criterion is minimized over the grid alone, so its minimum may lie beyond that end or between the end and its neighbor. Consider extending the grid past it and using a finer partition (more values), and compare the selected mediators and p-values across grids.
#>  term                  value   
#>  stat_hdmm             1.302   
#>  pval_hdmm             0.2538  
#>  stat_pe               865.5   
#>  j_pe                  864.2   
#>  pval_pe               < 0.0001
#>  total_indirect_effect -0.08396
#>  n_selected_mediators  3       
#>  df                    1       
#>  n_candidate_mediators 60      
#>  n_observations        200     
#> 
#> Outcome model: continuous (linear)
#> Selected mediators (3): 2, 3, 4
#> Tuning parameter (HBIC): lambda = 0.05 from 100 values in [0.05, 10] (the grid's lower end)
# HBIC chose 0.05, the smallest value of the default grid, so the fit
# warned (class poemed_grid_boundary). A grid that extends below it and
# is spaced more finely moves the choice inside the grid and keeps the
# same three mediators.
finer <- pe_mediation_linear(d$X, d$Y, d$M,
                             lambda_grid = seq(0.01, 10, length.out = 200))
attr(finer, "tuning")$lambda_selected
#> [1] 0.06020101
attr(finer, "selected_mediators")
#> [1] 2 3 4
```
