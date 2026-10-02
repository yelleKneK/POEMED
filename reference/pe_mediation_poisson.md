# Power-Enhanced Mediation Test for a Count Outcome

Tests the global null hypothesis of no active mediator for a count
outcome modeled with Poisson (log-link) regression, reporting both the
benchmark Wald test on the total indirect effect and the power-enhanced
(PE) test of Yu and Kelley (in press). This is the count-outcome worker
behind
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md).

## Usage

``` r
pe_mediation_poisson(
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

A tidy `data.frame` of class `poemed_tbl` with the same rows and
attributes as
[`pe_mediation_linear()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_linear.md),
except that the `"tuning"` list has no reduced-model entry.

## Details

The outcome follows a Poisson mediation model, \\Y \mid M, X, Z \sim
\mathrm{Poisson}(\exp(\alpha_m' M + \alpha_x' X + \alpha_z' Z))\\, with
a linear mediator model as in
[`pe_mediation_linear()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_linear.md).
Estimation, tuning, and the benchmark Wald test follow the
generalized-outcome framework of Guo et al. (2024): the mediator
coefficients are estimated by partial penalized likelihood with a SCAD
penalty (Fan and Li, 2001), leaving the exposures and confounders
unpenalized, with the tuning parameter chosen by the high-dimensional
BIC of Wang, Kim, and Li (2013). The penalized Poisson fit is a glmnet
lasso fit followed by one weighted lasso step on a quadratic
approximation of the Poisson likelihood at the lasso estimate, with the
SCAD derivatives at the lasso coefficients as weights: the one-step
local linear approximation of Zou and Li (2008), started from the lasso
as in Fan, Xue, and Zou (2014). The power-enhanced test of Yu and Kelley
(in press), built on the power enhancement principle of Fan, Liao, and
Yao (2015), adds the component \\J_m\\ (see
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)).
The total indirect effect \\\beta = \Gamma_x \alpha_m\\ is on the
log-mean (link) scale, not the count scale: with `scale = TRUE` it is
the change in the log mean count per standard deviation of the exposure.

## How to Cite

If you use POEMED in published work, please cite Yu and Kelley (in
press), the article that introduces its methods and the package.
`citation("POEMED")` gives the full reference and a BibTeX entry.

## References

Yu, X., & Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

Guo, X., Li, R., Liu, J., & Zeng, M. (2024). Estimations and tests for
generalized mediation models with high-dimensional potential mediators.
*Journal of Business & Economic Statistics, 42*(1), 243–256.
[doi:10.1080/07350015.2023.2174548](https://doi.org/10.1080/07350015.2023.2174548)

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

[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md),
[`pe_mediation_linear()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_linear.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_logistic.md).

Other mediation tests:
[`pe_mediate()`](https://yelleknek.github.io/POEMED/reference/pe_mediate.md),
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md),
[`pe_mediation_linear()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_linear.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_logistic.md),
[`pe_mediators()`](https://yelleknek.github.io/POEMED/reference/pe_mediators.md),
[`pe_selection()`](https://yelleknek.github.io/POEMED/reference/pe_selection.md),
[`summary.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/summary.poemed_tbl.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
                             outcome = "count", c1 = 0.5, c2 = 0.4)
pe_mediation_poisson(d$X, d$Y, d$M)
#>  term                  value   
#>  stat_hdmm             0.001175
#>  pval_hdmm             0.9727  
#>  stat_pe               545.4   
#>  j_pe                  545.4   
#>  pval_pe               < 0.0001
#>  total_indirect_effect 0.001457
#>  n_selected_mediators  2       
#>  df                    1       
#>  n_candidate_mediators 60      
#>  n_observations        200     
#> 
#> Outcome model: count (Poisson)
#> Selected mediators (2): 4, 5
#> Tuning parameter (HBIC): lambda = 0.251 from 100 values in [0.05, 10]
```
