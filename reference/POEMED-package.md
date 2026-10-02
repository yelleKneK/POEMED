# Power Enhancement for High-Dimensional Mediation Analysis (POEMED)

Powerful global tests for high-dimensional mediation analysis, for
continuous, binary, and count outcomes, implementing the power-enhanced
tests of Yu and Kelley (in press). The package answers the question “is
there any active mediator among a large set of candidate mediators?”
and, when the answer is yes, reports which individual mediators the
screen selects. The name comes from POwer-Enhanced MEDiation and is
pronounced “POE-med”.

## The Problem POEMED Solves

A mediation analysis asks how an exposure \\X\\ affects an outcome \\Y\\
through intermediate variables (mediators) \\M\\. With many candidate
mediators, the standard approach tests the *total* indirect effect
\\\beta = \Gamma_x \alpha_m\\, the sum of every mediator's individual
indirect effect, with the partially penalized Wald test of Guo et al.
(2022, 2023, 2024). That test has a blind spot: when some mediators
carry positive indirect effects and others carry negative ones, the
individual effects can cancel and the total indirect effect is zero even
though active mediators plainly exist. A test built on the total
indirect effect is then powerless, by construction, against exactly the
heterogeneous mediation that is common in practice.

The power enhancement (PE) tests of Yu and Kelley (in press), which
follow the power enhancement principle of Fan, Liao, and Yao (2015), add
a component \\J_m\\ that accumulates the marginal signal from each
individual mediator. Because \\J_m\\ is a sum of magnitudes,
opposite-signed indirect effects reinforce rather than cancel, so the PE
test detects active mediators whether their effects are homogeneous,
heterogeneous, or exactly contrasting. Under the global null of no
mediation the PE test has the same chi-square limit as the original
test, so its Type I error rate is controlled asymptotically; in finite
samples it can reject somewhat more often than the original test (see
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)).

## Main Functions

- [`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md):

  The front end. Give it the exposure, outcome, mediators, and optional
  confounders, name the outcome type (continuous, binary, or count), and
  it returns a tidy table with the benchmark Wald test, the
  power-enhanced test, the estimated total indirect effect, and the set
  of mediators the screen selected.

- [`pe_mediation_linear()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_linear.md),
  [`pe_mediation_logistic()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_logistic.md),
  [`pe_mediation_poisson()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_poisson.md):

  The outcome-specific workers that
  [`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
  dispatches to. Call them directly when the outcome type is known.

- [`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md):

  Generates mediation data under the homogeneous and contrasting
  patterns studied in the article, for any of the three outcome types.

- [`simulate_guo_mediation()`](https://yelleknek.github.io/POEMED/reference/simulate_guo_mediation.md)
  and
  [guo_calibration](https://yelleknek.github.io/POEMED/reference/guo_calibration.md):

  The article's real-data-motivated heterogeneous setting (Section S.4.4
  of its supplement), a continuous-outcome design with eleven active
  mediators of mixed sign, and the constants that calibrate it.

- [`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md):

  Runs the article's Monte Carlo size and power study: it sweeps a
  signal-strength grid and returns the empirical rejection rate of the
  benchmark and power-enhanced tests at each point.

- [`pe_mediate()`](https://yelleknek.github.io/POEMED/reference/pe_mediate.md):

  A data-frame / formula front end to
  [`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
  for users who prefer naming columns to matrices.

- [`pe_mediators()`](https://yelleknek.github.io/POEMED/reference/pe_mediators.md)
  and
  [`pe_selection()`](https://yelleknek.github.io/POEMED/reference/pe_selection.md):

  The per-mediator screen detail behind a fit, and the selected set
  under each multiplicity method.

- [`pe_simulation_study()`](https://yelleknek.github.io/POEMED/reference/pe_simulation_study.md)
  and
  [`ss_power_pe_mediation()`](https://yelleknek.github.io/POEMED/reference/ss_power_pe_mediation.md):

  Run the size and power study across patterns, and plan the sample size
  for a target power.

- [`pe_identification_study()`](https://yelleknek.github.io/POEMED/reference/pe_identification_study.md):

  Runs the article's Monte Carlo study of mediator identification: the
  familywise error rate, false discovery rate, precision, and recall of
  the selected set under each multiplicity method.

- [`plot.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/plot.poemed_tbl.md),
  [`tidy.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/poemed_broom.md),
  [`glance.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/poemed_broom.md):

  A base-graphics plot and broom verbs for POEMED results.

- [`mediation_ar1_cov()`](https://yelleknek.github.io/POEMED/reference/mediation_ar1_cov.md):

  The autoregressive mediator covariance used throughout the
  simulations, exported as a reusable utility.

- [WHO_health_mediation](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md),
  [`WHO_mediation_design()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_design.md),
  and
  [`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md):

  The article's empirical data (economic growth, 57 health-expenditure
  mediators, five health outcomes for 91 WHO members), a helper that
  assembles an analysis design from it, and a function that runs the
  article's subgroup analyses (global, by region, by income group, and
  by region and income group) with the preprocessing and the tuning grid
  of the analysis behind the article's tables (standardized exposure and
  mediators, centered outcome, unscaled region, income, and year
  confounders, and the grid `seq(0.1, 2, length.out = 100)`; the
  article's text says that covariates are standardized). Its country
  counts match the article's tables in every cell, except that the
  function reports `NA` for the cells the article lists with 0
  countries. Under its default grid it reproduces the article's Tables
  1, 2, and S.8 to S.12 in every cell it fits: every selected set but
  one BH set, and the benchmark p-values to four decimals in all but a
  few cells. Six region-by-income cells are reported as `NA` instead;
  the reproduction vignette,
  [`vignette("reproducing-the-article")`](https://yelleknek.github.io/POEMED/articles/reproducing-the-article.md),
  lists those cells and the ones that differ, and why.

## Reading the Output

Every estimation and testing function returns a tidy `data.frame` with a
`term` column and a numeric `value` column, carrying the `poemed_tbl`
class so it prints with whole numbers shown without a decimal part,
other quantities to a few significant figures, and p-values to four
decimal places. The stored numbers keep full precision; only the display
rounds. See
[poemed_tbl](https://yelleknek.github.io/POEMED/reference/poemed_tbl.md)
for the details.

## How to Cite

If you use POEMED in published work, please cite Yu and Kelley (in
press), the article that introduces its methods and the package.
`citation("POEMED")` gives the full reference and a BibTeX entry.

## References

Yu, X., & Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

Fan, J., Liao, Y., & Yao, J. (2015). Power enhancement in
high-dimensional cross-sectional tests. *Econometrica, 83*(4),
1497–1541. [doi:10.3982/ECTA12749](https://doi.org/10.3982/ECTA12749)

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

Guo, X., Li, R., Liu, J., & Zeng, M. (2024). Estimations and tests for
generalized mediation models with high-dimensional potential mediators.
*Journal of Business & Economic Statistics, 42*(1), 243–256.
[doi:10.1080/07350015.2023.2174548](https://doi.org/10.1080/07350015.2023.2174548)

## See also

Useful links:

- <https://github.com/yelleKneK/POEMED>

- <https://yelleknek.github.io/POEMED/>

- Report bugs at <https://github.com/yelleKneK/POEMED/issues>

## Author

Xiufan Yu and Ken Kelley
