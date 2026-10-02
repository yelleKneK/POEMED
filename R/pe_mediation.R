#' Power-Enhanced Test for High-Dimensional Mediation
#'
#' Tests whether any mediator is active among a large, possibly
#' intercorrelated set of candidate mediators, for a continuous, binary,
#' or count outcome, with the power-enhanced (PE) test of Yu and Kelley
#' (in press). This is the main entry point of \pkg{POEMED}: name the
#' outcome type and it dispatches to the appropriate model, returning one
#' tidy table that reports the benchmark Wald test on the total indirect
#' effect side by side with the PE test, plus the set of individual
#' mediators the screen selected.
#'
#' @details
#' Standard high-dimensional mediation tests target the \emph{total}
#' indirect effect \eqn{\beta = \Gamma_x \alpha_m}{beta = Gamma_x alpha_m},
#' the sum of every mediator's individual indirect effect. When some indirect
#' effects are positive and others negative they can cancel, making
#' \eqn{\beta = 0}{beta = 0} even though active mediators exist, and a test
#' built on \eqn{\beta}{beta} is then powerless. The benchmark reported here
#' is such a test, the Wald test \eqn{S_n}{S_n} on \eqn{\beta}{beta}. For a
#' continuous outcome it is the test of Guo et al. (2022), which extends the
#' partially penalized Wald test of Guo et al. (2023) to observed
#' confounders. For binary and count outcomes it is the test of Guo et al.
#' (2024).
#'
#' The PE test of Yu and Kelley (in press) follows the power enhancement
#' principle of Fan, Liao, and Yao (2015). It adds to \eqn{S_n}{S_n} the
#' component (equation 2.9 of the article)
#' \deqn{J_m = \sqrt{p} \sum_{i=1}^{q} \sum_{j \in \hat{S}}
#'   \left| \frac{\hat\alpha_{m,j}}{\hat\sigma_{m,j}} \right|
#'   \left| \frac{\hat\Gamma_{x,i,j}}{\hat\sigma_{\Gamma,i,j}} \right|
#'   \mathbf{1}\!\left\{ \max(p_{1,j}, p_{2,i,j}) <
#'   \frac{\alpha_{\mathrm{lvl}}}{s q \log\log n} \right\},}{
#'   J_m = sqrt(p) sum_(i = 1, ..., q) sum_(j in hat(S))
#'   |hat(alpha)_m[j] / hat(sigma)_m[j]|
#'   |hat(Gamma)_x[i, j] / hat(sigma)_Gamma[i, j]|
#'   1{max(p_1[j], p_2[i, j]) < alpha_lvl / (s q log(log(n)))},}
#' a sum over the \eqn{s}{s} selected mediators \eqn{\hat S}{hat(S)} of the
#' product of the (standardized) mediator-on-outcome and exposure-on-mediator
#' statistics, kept only for pairs where both paths are individually
#' significant; \eqn{\alpha_{\mathrm{lvl}}}{alpha_lvl} is `error_level`. That
#' indicator is the Bonferroni screen used by `method = "Bonferroni"`.
#' Under `method = "BH"` or `"BY"` the indicator is instead
#' \eqn{\mathbf{1}\{\tilde p_{i,j} < \alpha_{\mathrm{lvl}} / \log\log n\}}{
#' 1{tilde(p)[i, j] < alpha_lvl / log(log(n))}},
#' where \eqn{\tilde p_{i,j}}{tilde(p)[i, j]} is the Benjamini and Hochberg
#' (1995) or Benjamini and Yekutieli (2001) adjusted value of
#' \eqn{\max(p_{1,j}, p_{2,i,j})}{max(p_1[j], p_2[i, j])} over the
#' \eqn{s q}{s q} pairs (equations S.22 to S.25 of the article's supplement).
#' Because \eqn{J_m}{J_m} accumulates magnitudes, opposite-signed indirect
#' effects reinforce rather than cancel, so the test
#' \eqn{M_{PE} = S_n + J_m}{M_PE = S_n + J_m} stays powerful under
#' heterogeneous and contrasting mediation. The indicator also selects
#' individual mediators, with asymptotic familywise error rate control
#' under `method = "Bonferroni"` or false discovery rate control under
#' `"BH"` / `"BY"`. In this documentation an \emph{active} mediator is one
#' whose effect is truly nonzero (the article's usage, as in the global
#' null of no active mediator), and a \emph{selected} mediator is one the
#' screen picks as statistically significant: the selected set is the
#' method's estimate of the active set.
#'
#' Both statistics are referred to a chi-square distribution with
#' \eqn{q}{q} degrees of freedom. Under the global null of no active
#' mediator, \eqn{J_m}{J_m} is zero with probability tending to one, so
#' \eqn{M_{PE}}{M_PE} has the same chi-square limit as \eqn{S_n}{S_n} and the
#' PE test is asymptotically valid for any fixed `error_level`, by Theorem 1
#' of Yu and Kelley (in press). That guarantee is asymptotic. In finite
#' samples, under a null in which some mediators affect the outcome but the
#' exposure affects none of them, the screen passes one of those mediators
#' with probability of about `error_level / log(log(n))`, which is 0.030 at
#' \eqn{n = 200}{n = 200} with `error_level = 0.05`, and \eqn{J_m}{J_m} is
#' then large enough to reject. The PE test's size can therefore exceed the
#' benchmark's by up to about that amount. For example,
#' `pe_power_curve(200, 100, pattern = "homogeneous", c1_grid = 0,
#' n_rep = 2000, seed = 113)` gave rejection rates at the 0.05 level of
#' 0.0475 for the benchmark and 0.0620 for the PE test over 2,000 null
#' samples, and 0.0510 for the PE test with `error_level = 0.01` (Monte
#' Carlo standard errors about 0.005). A smaller
#' `error_level` narrows the gap, at some cost in identifying active
#' mediators.
#'
#' @param X Numeric matrix of exposures with \eqn{n}{n} rows and \eqn{q}{q}
#'   columns (\eqn{q}{q} is usually 1). A numeric vector is treated as a
#'   single exposure.
#' @param Y Numeric outcome vector of length \eqn{n}{n}. Continuous for
#'   `outcome = "continuous"`, coded 0/1 for `"binary"`, and nonnegative
#'   integer counts for `"count"`.
#' @param M Numeric matrix of candidate mediators with \eqn{n}{n} rows and
#'   \eqn{p}{p} columns; \eqn{p}{p} may exceed \eqn{n}{n}. Column names, when
#'   present, label the selected mediators in the print footer and the
#'   mediator table; without them mediators are reported by column
#'   position.
#' @param Z Optional numeric matrix of confounders with \eqn{n}{n} rows.
#'   Default `NULL` (no confounders).
#' @param outcome Outcome type: `"continuous"` (linear model),
#'   `"binary"` (logistic model), or `"count"` (Poisson log-link model).
#' @param method Multiplicity adjustment used to screen individual
#'   mediators in the PE component: `"Bonferroni"` (familywise error
#'   rate control, the default), `"BH"`, or `"BY"` (false discovery rate
#'   control by the procedures of Benjamini and Hochberg, 1995, and of
#'   Benjamini and Yekutieli, 2001; `"BH"` assumes independence or
#'   positive dependence among mediators, `"BY"` is valid under arbitrary
#'   dependence). See Details for the screening threshold each one uses.
#' @param scale Logical; if `TRUE` (default) the exposures, mediators, and
#'   confounders are each standardized to mean 0 and standard deviation 1
#'   before fitting, as in the article, and a continuous outcome is
#'   centered but not rescaled. The total indirect effect is then the
#'   change in the outcome per standard deviation of each exposure: in the
#'   units of `Y` for a continuous outcome, in log odds for a binary
#'   outcome, and in log mean count for a count outcome. With `FALSE` the
#'   inputs are used exactly as supplied, and the total indirect effect is
#'   per unit of each exposure. The continuous-outcome fit has no
#'   intercept, and the default tuning grid was chosen for standardized
#'   columns, so `FALSE` is meant for inputs that are already standardized (and,
#'   for a continuous outcome, a `Y` that is already centered); on raw
#'   inputs it fits a different model and can even change the sign of the
#'   estimate. Because a continuous `Y` is not rescaled, the tuning grid is
#'   not expressed in the outcome's units, so the selected mediators, both
#'   p-values, and the verdict depend on the units of `Y`: the same data
#'   with `Y` multiplied by 10 can select a different set of mediators.
#' @param error_level Target error rate for the mediator-screening step,
#'   in \eqn{(0, 1)}{(0, 1)}. Interpreted as the familywise error rate when
#'   `method = "Bonferroni"` and as the false discovery rate when
#'   `method = "BH"` or `"BY"`. Default 0.05. This is distinct from the
#'   significance level used to test the global null (which is the user's
#'   choice when reading `pval_pe`). The PE test is asymptotically valid
#'   for any fixed value in \eqn{(0, 1)}{(0, 1)}, but in finite samples a larger
#'   value lets the screen pass more spurious mediators under the null and
#'   raises the PE test's size further above the benchmark's (see
#'   Details).
#' @param report_all_methods Logical; if `TRUE`, the selected-mediator set is
#'   computed for all three multiplicity methods (Bonferroni, BH, BY) off the
#'   same fit and recorded for comparison, retrievable with [pe_selection()].
#'   The `method` argument still drives the primary result. Default `FALSE`.
#' @param drop_constant Logical; how to handle constant (zero-variance)
#'   mediator columns, which cannot be standardized or carry signal. If
#'   `FALSE` (the default) they are an error; if `TRUE` they are dropped
#'   with a warning and the remaining mediators are renumbered back to their
#'   original column positions in the reported selected set. Useful for small
#'   subgroups where some mediators happen to be constant.
#' @param lambda_grid Numeric vector of candidate tuning parameters for the
#'   SCAD penalty (Fan and Li, 2001) in the penalized mediator fit, which
#'   for a binary outcome is adaptively rescaled (see
#'   [pe_mediation_logistic()]). The fit is repeated at each value and the
#'   one minimizing the high-dimensional BIC (HBIC) of Wang, Kim, and Li
#'   (2013) is kept: whatever the grid, the value selected is the best of
#'   the values offered, since HBIC is minimized over the grid and not over
#'   every positive lambda. The default, `seq(0.05, 10, length.out = 100)`
#'   for every outcome family, is written out in the signature so that it
#'   can be edited in place. Every value is fit (the search is not adaptive),
#'   so the cost of a fit grows with the length of the grid. The grid
#'   can be tuned to the situation. Fewer values, or a range narrowed to
#'   where earlier fits put the HBIC minimum, make each fit faster, which
#'   matters when many fits are run (a simulation study, for example); a
#'   wider or finer grid is worth trying when the selected value is the
#'   smallest or largest value of the grid, in which case a non-empty fit
#'   warns (class `poemed_grid_boundary`) and the print footer, which
#'   reports the full model's choice, says so. For a continuous outcome
#'   the benchmark Wald test also needs a penalized fit of the reduced
#'   model, the outcome on the mediators and confounders without the
#'   exposure; that fit searches the same grid, its choice is
#'   `lambda_selected_reduced` in the `"tuning"` attribute, and the warning
#'   covers it as well, naming whether the full model, the reduced model,
#'   or both sat at the end. For a
#'   continuous outcome the grid is fixed while `Y` is centered but not
#'   rescaled, so the value HBIC chooses, and with it the selected
#'   mediators, depends on the units of `Y` (see `scale`). The grid
#'   searched and the value chosen are reported in the `"tuning"`
#'   attribute and the print footer.
#'
#' @return A tidy `data.frame` of class `poemed_tbl`. See
#'   [pe_mediation_linear()] for the row schema and the attributes. The
#'   selected mediators (the candidates the screen picked; see Details for
#'   the distinction from active mediators) are in the
#'   `"selected_mediators"` attribute and the print footer; the
#'   tuning parameter HBIC chose is in the `"tuning"` attribute and the
#'   footer. The total indirect effect is on the scale described under
#'   `scale`.
#'
#' @section When to Use POEMED:
#' POEMED is for the global question, \dQuote{is there any active mediator
#' among many candidates?}, followed by the selection of individual
#' mediators. It assumes a linear or generalized linear outcome model, a sparse
#' set of active mediators, standardized inputs, and complete data. Some
#' practical limits, seen on simulated data:
#' \itemize{
#'   \item Sample size. The penalized selection needs \eqn{n}{n} large
#'     relative to \eqn{\log p}{log(p)} and to the signal. With
#'     \eqn{n = 100}{n = 100} and \eqn{p = 2000}{p = 2000} under the contrasting
#'     pattern (continuous outcome, \eqn{c_1 = 1}{c1 = 1}, default grid, seeds 1
#'     to 10) the selection kept at most one mediator: one in 4 of the 10
#'     samples and none in the other 6. The power-enhanced test could then
#'     draw on that one mediator only, not on the contrasting structure; it
#'     rejected in 4 of the 10 samples, and the benchmark in 2.
#'   \item Strongly correlated mediators. Correlation among the mediators
#'     weakens both the selection and the screen. Under the contrasting
#'     pattern (continuous outcome, \eqn{n = 200}{n = 200},
#'     \eqn{p = 60}{p = 60}, \eqn{c_1 = 1}{c1 = 1}, default grid, seeds 1 to 10)
#'     the screen selected at least one mediator in 6 of 10 samples with an
#'     autoregressive correlation of 0.9 among neighboring mediators, against
#'     10 of 10 at 0.5, and the selected set varied from sample to sample
#'     (one sample named a null mediator beside an active one at both
#'     correlations).
#'   \item Selected mediators in finite samples. The screen's control of
#'     its error rate is asymptotic. Under the contrasting pattern with a
#'     strong signal (continuous outcome, \eqn{c_1 = 1}{c1 = 1}, default
#'     grid, Bonferroni screen, `error_level = 0.05`) the selected set held
#'     an inactive mediator in 0.14 of 400 samples at \eqn{n = 150}{n = 150},
#'     \eqn{p = 60}{p = 60} and in 0.11 of 200 samples at
#'     \eqn{n = 300}{n = 300}, \eqn{p = 500}{p = 500}; the rate was 0.05
#'     under the homogeneous pattern at the smaller design and 0.00 with no
#'     signal. It depends on the tuning grid: 20 values from 0.2 to 0.39 gave
#'     0.00 in the same runs, with lower recall (0.13 against 0.48 at the
#'     smaller design, 0.25 against 0.47 at the larger).
#'   \item Rare binary outcomes. With about 2 to 3 percent events at
#'     \eqn{n = 300}{n = 300} the penalized logistic fit selects nothing at any
#'     grid value; the fit is empty and the function says so.
#'   \item Overdispersed counts. The count model is Poisson. In null
#'     simulations with negative binomial counts (size parameter 1) at
#'     \eqn{n = 200}{n = 200}, \eqn{p = 100}{p = 100}, 300 replications on a
#'     grid of 20 values from 0.05 to 1, its size was not inflated, but about
#'     a tenth of the fits were empty. A count
#'     fit can also take several times as long as a continuous fit of the same
#'     design.
#'   \item Missing values are not handled; supply complete cases.
#' }
#' For a single mediator, or a few mediators fit as one structural model
#' (indirect effects with confidence intervals, likelihood ratio tests of
#' arbitrary indirect effects, moderated mediation), the \pkg{DMAR}
#' package is the tool; POEMED's contribution is the high-dimensional global
#' test. The other high-dimensional methods the article compares against
#' are HILMA (Zhou, Wang, & Zhao, 2020), GlobalTest (Djordjilovic et al.,
#' 2019), HDMT (Dai, Stanford, & LeBlanc, 2022), and DACT (Liu et al.,
#' 2022). POEMED does not depend on them. HDMT is on CRAN, and GlobalTest
#' is the Bioconductor package \pkg{globaltest}. The CRAN package
#' \pkg{freebird}, which implemented HILMA, was archived in July 2026, and
#' the CRAN package named \pkg{DACT} is an unrelated clinical-trials
#' package, not the method of Liu et al. The vignette
#' `poemed-vs-competitors` shows how to call HILMA and GlobalTest beside a
#' POEMED fit.
#'
#' @references
#' Yu, X., & Kelley, K. (in press). Power Enhancement in
#' High-Dimensional Heterogeneous Mediation Analysis. \emph{Journal of the
#' American Statistical Association}.
#'
#' Fan, J., Liao, Y., & Yao, J. (2015). Power enhancement in
#' high-dimensional cross-sectional tests. \emph{Econometrica, 83}(4),
#' 1497--1541. \doi{10.3982/ECTA12749}
#'
#' Guo, X., Li, R., Liu, J., & Zeng, M. (2022). High-dimensional
#' mediation analysis for selecting DNA methylation loci mediating
#' childhood trauma and cortisol stress reactivity. \emph{Journal of the
#' American Statistical Association, 117}(539), 1110--1121.
#' \doi{10.1080/01621459.2022.2053136}
#'
#' Guo, X., Li, R., Liu, J., & Zeng, M. (2023). Statistical inference for
#' linear mediation models with high-dimensional mediators and application
#' to studying stock reaction to COVID-19 pandemic. \emph{Journal of
#' Econometrics, 235}(1), 166--179. \doi{10.1016/j.jeconom.2022.03.001}
#'
#' Guo, X., Li, R., Liu, J., & Zeng, M. (2024). Estimations and tests
#' for generalized mediation models with high-dimensional potential
#' mediators. \emph{Journal of Business & Economic Statistics, 42}(1),
#' 243--256. \doi{10.1080/07350015.2023.2174548}
#'
#' Fan, J., & Li, R. (2001). Variable selection via nonconcave penalized
#' likelihood and its oracle properties. \emph{Journal of the American
#' Statistical Association, 96}(456), 1348--1360.
#' \doi{10.1198/016214501753382273}
#'
#' Wang, L., Kim, Y., & Li, R. (2013). Calibrating nonconvex penalized
#' regression in ultra-high dimension. \emph{The Annals of Statistics,
#' 41}(5), 2505--2536. \doi{10.1214/13-AOS1159}
#'
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
#' rate: A practical and powerful approach to multiple testing.
#' \emph{Journal of the Royal Statistical Society, Series B, 57}(1),
#' 289--300. \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Benjamini, Y., & Yekutieli, D. (2001). The control of the false
#' discovery rate in multiple testing under dependency. \emph{The Annals
#' of Statistics, 29}(4), 1165--1188. \doi{10.1214/aos/1013699998}
#'
#' Zhou, R. R., Wang, L., & Zhao, S. D. (2020). Estimation and inference
#' for the indirect effect in high-dimensional linear mediation models.
#' \emph{Biometrika, 107}(3), 573--589. \doi{10.1093/biomet/asaa016}
#'
#' Djordjilovic, V., Page, C. M., Gran, J. M., Nost, T. H., Sandanger,
#' T. M., Veierod, M. B., & Thoresen, M. (2019). Global test for
#' high-dimensional mediation: Testing groups of potential mediators.
#' \emph{Statistics in Medicine, 38}(18), 3346--3360.
#' \doi{10.1002/sim.8199}
#'
#' Dai, J. Y., Stanford, J. L., & LeBlanc, M. (2022). A multiple-testing
#' procedure for high-dimensional mediation hypotheses. \emph{Journal of
#' the American Statistical Association, 117}(537), 198--213.
#' \doi{10.1080/01621459.2020.1765785}
#'
#' Liu, Z., Shen, J., Barfield, R., Schwartz, J., Baccarelli, A. A., &
#' Lin, X. (2022). Large-scale hypothesis testing for causal mediation
#' effects with applications in genome-wide epigenetic studies.
#' \emph{Journal of the American Statistical Association, 117}(537),
#' 67--81. \doi{10.1080/01621459.2021.1914634}
#'
#' @inheritSection POEMED-package How to Cite
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso The outcome-specific workers [pe_mediation_linear()],
#'   [pe_mediation_logistic()], [pe_mediation_poisson()];
#'   [simulate_mediation_data()] to generate data; [pe_power_curve()] to
#'   run the designs of the article's size and power studies.
#'
#' @family mediation tests
#'
#' @examples
#' set.seed(113)
#' # Contrasting mediation: active mediators whose effects cancel, so the
#' # total indirect effect is zero. Compare the two p-values: the benchmark
#' # Wald test (pval_hdmm) is large while the PE test (pval_pe) is tiny.
#' d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
#'                              outcome = "continuous", c1 = 1)
#' pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
#' # HBIC chose 0.05, the smallest value of the default grid, so the fit
#' # warned (class poemed_grid_boundary). A grid that extends below it and
#' # is spaced more finely moves the choice inside the grid and keeps the
#' # same three mediators.
#' finer <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
#'                       lambda_grid = seq(0.01, 10, length.out = 200))
#' attr(finer, "tuning")$lambda_selected
#' attr(finer, "selected_mediators")
#'
#' @export
pe_mediation <- function(X, Y, M, Z = NULL,
                         outcome = c("continuous", "binary", "count"),
                         method = c("Bonferroni", "BH", "BY"),
                         scale = TRUE, error_level = 0.05,
                         report_all_methods = FALSE,
                         drop_constant = FALSE,
                         lambda_grid = seq(0.05, 10, length.out = 100)) {
  outcome <- .poemed_match_arg(outcome)
  method <- .poemed_match_arg(method)
  worker <- switch(outcome,
                   continuous = pe_mediation_linear,
                   binary     = pe_mediation_logistic,
                   count      = pe_mediation_poisson)
  worker(X = X, Y = Y, M = M, Z = Z, method = method, scale = scale,
         error_level = error_level,
         report_all_methods = report_all_methods, drop_constant = drop_constant,
         lambda_grid = lambda_grid)
}
