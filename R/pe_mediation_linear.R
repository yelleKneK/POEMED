#' Power-Enhanced Mediation Test for a Continuous Outcome
#'
#' Tests the global null hypothesis of no active mediator among a
#' high-dimensional set of candidate mediators, for a continuous
#' (linear-model) outcome, and reports both the benchmark Wald test on
#' the total indirect effect and the power-enhanced (PE) test of Yu and
#' Kelley (in press), which remains powerful when individual mediation
#' effects are heterogeneous or contrasting. This is the continuous-outcome
#' worker behind [pe_mediation()]; call it directly when the outcome is
#' continuous.
#'
#' @details
#' The model is the linear mediation pair
#' \deqn{Y = \alpha_m' M + \alpha_x' X + \alpha_z' Z + \varepsilon_y,
#'       \qquad M = \Gamma_x' X + \Gamma_z' Z + \varepsilon_m,}{
#'   Y = alpha_m' M + alpha_x' X + alpha_z' Z + eps_y,
#'   M = Gamma_x' X + Gamma_z' Z + eps_m,}
#' with total indirect effect
#' \eqn{\beta = \Gamma_x \alpha_m}{beta = Gamma_x alpha_m}.
#' Following Guo et al. (2023), the mediator coefficients
#' \eqn{\alpha_m}{alpha_m} are estimated by partial penalized least squares
#' with a SCAD penalty (Fan and Li, 2001), leaving the exposures and
#' confounders unpenalized. The penalized fit is computed by the one-step
#' local linear approximation of Zou and Li (2008) started from a lasso fit,
#' as in Fan, Xue, and Zou (2014): a lasso fit, then a weighted lasso fit
#' whose weights are the SCAD derivatives at the lasso coefficients. The
#' tuning parameter is chosen by the high-dimensional BIC of Wang, Kim, and
#' Li (2013). The fit has no intercept, because with `scale = TRUE` the
#' outcome is centered and the other inputs standardized first (see `scale`
#' in [pe_mediation()]).
#'
#' The benchmark statistic is the Wald statistic
#' \eqn{S_n = n \hat\beta' \hat\Sigma_\beta^{-1} \hat\beta}{
#' S_n = n hat(beta)' hat(Sigma)_beta^(-1) hat(beta)}
#' of Guo et al. (2022), which extends the partially penalized Wald test of
#' Guo et al. (2023) to observed confounders. The power-enhanced statistic of
#' Yu and Kelley (in press), built on the power enhancement principle of Fan,
#' Liao, and Yao (2015), adds the component \eqn{J_m}{J_m} that accumulates
#' the marginal signal from each selected mediator (see [pe_mediation()] for
#' the rationale, the formula, the screening thresholds, and the test's
#' finite-sample size). Both are referred to a chi-square distribution with
#' \eqn{q}{q} (the number of exposures) degrees of freedom, their common
#' limit under the global null.
#'
#' @inheritParams pe_mediation
#'
#' @return A tidy `data.frame` of class `poemed_tbl` with rows
#'   `stat_hdmm` and `pval_hdmm` (the benchmark Wald test),
#'   `stat_pe`, `j_pe`, and `pval_pe` (the power-enhanced test and its
#'   PE component), `total_indirect_effect` (the estimated total indirect
#'   effect; one row per exposure when `q > 1`, suffixed with the
#'   exposure's column name when `X` has column names and with `_1`, `_2`,
#'   ... otherwise), `n_selected_mediators` (the number of mediators the
#'   screen selected), `df` (the chi-square degrees of freedom,
#'   equal to `q`), `n_candidate_mediators`, and `n_observations`. With
#'   `scale = TRUE` the total indirect effect is in the units of `Y` per
#'   standard deviation of the exposure; with `scale = FALSE` it is in the
#'   units of `Y` per unit of the exposure as supplied. No confidence
#'   interval is reported: the article reports none, and the Wald interval
#'   built on the penalized estimate under-covers in finite samples.
#'
#'   Attributes: `"selected_mediators"` (the column positions in `M` that
#'   the screen selected; the print footer shows their column
#'   names when `M` has them), `"mediator_names"` and `"exposure_names"` (the column
#'   names of `M` and `X`, or `NULL` when unnamed), `"mediator_table"`
#'   (the per-mediator screening statistics for the penalized-selected
#'   mediators, see [pe_mediators()]; its `t_exposure` is the signed
#'   exposure-on-mediator statistic when `q = 1` and the largest absolute
#'   one across exposures when `q > 1`, and its `screen_p` is the smallest
#'   screening p-value across exposures), `"method"`, `"error_level"`,
#'   `"outcome"`, `"p_terms"` (the names of the rows that
#'   hold p-values, which the print method formats as p-values),
#'   `"tuning"` (a list with the grid searched, the HBIC-selected
#'   `lambda_selected`, and `at_lower_end` and `at_upper_end`, `TRUE` when
#'   the selection sat at the grid's smallest or largest value; for a
#'   continuous outcome also `lambda_selected_reduced`, the choice for the
#'   reduced model of the benchmark Wald test, searched on the same grid),
#'   `"empty_fit"`, and `"pe_by_method"` (a data frame with one row per
#'   multiplicity method screened, giving that method's `stat_pe`, `j_pe`,
#'   `pval_pe`, and number of selected mediators; it has one row unless
#'   `report_all_methods = TRUE`). With `report_all_methods = TRUE` there is
#'   also `"selection_by_method"`, the selected set under each method (see
#'   [pe_selection()]).
#'
#'   When the penalized fit that HBIC selects contains no mediator there
#'   is nothing to test: the function warns (a condition of class
#'   `poemed_empty_fit`), sets `"empty_fit"` to `TRUE`, and returns the same
#'   rows with both statistics 0, both p-values 1 (POEMED's reporting
#'   convention, under which an empty selection is a non-rejection; the
#'   article states only that \eqn{J_m = 0}{J_m = 0} when the selection is
#'   empty), and a total indirect effect of 0. A
#'   p-value of 1 from such a fit is not evidence for the null. An empty fit is
#'   expected when no mediator is active. A non-empty fit warns (class
#'   `poemed_grid_boundary`) when its HBIC-selected lambda is the smallest or
#'   largest value of the grid, for the full model or, for a continuous
#'   outcome, for the reduced model of the benchmark Wald test, and the
#'   warning names the fit that sat there: the criterion is minimized over
#'   the grid alone, so extend the grid past that end and space its values
#'   more finely.
#'
#' @references
#' Yu, X., & Kelley, K. (in press). Power Enhancement in
#' High-Dimensional Heterogeneous Mediation Analysis. \emph{Journal of the
#' American Statistical Association}.
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
#' Fan, J., Liao, Y., & Yao, J. (2015). Power enhancement in
#' high-dimensional cross-sectional tests. \emph{Econometrica, 83}(4),
#' 1497--1541. \doi{10.3982/ECTA12749}
#'
#' Fan, J., & Li, R. (2001). Variable selection via nonconcave penalized
#' likelihood and its oracle properties. \emph{Journal of the American
#' Statistical Association, 96}(456), 1348--1360.
#' \doi{10.1198/016214501753382273}
#'
#' Zou, H., & Li, R. (2008). One-step sparse estimates in nonconcave
#' penalized likelihood models. \emph{The Annals of Statistics, 36}(4),
#' 1509--1533. \doi{10.1214/009053607000000802}
#'
#' Fan, J., Xue, L., & Zou, H. (2014). Strong oracle optimality of folded
#' concave penalized estimation. \emph{The Annals of Statistics, 42}(3),
#' 819--849. \doi{10.1214/13-AOS1198}
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
#' @inheritSection POEMED-package How to Cite
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()] for the front end, [pe_mediation_logistic()]
#'   and [pe_mediation_poisson()] for binary and count outcomes, and
#'   [simulate_mediation_data()] to generate data for trying the method.
#'
#' @family mediation tests
#'
#' @examples
#' set.seed(113)
#' # A contrasting setting: four active mediators whose indirect effects
#' # cancel, so the total indirect effect is zero. The Wald test does not
#' # reject here; the PE test rejects and selects three of the four
#' # active mediators.
#' d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
#'                              outcome = "continuous", c1 = 1)
#' pe_mediation_linear(d$X, d$Y, d$M)
#' # HBIC chose 0.05, the smallest value of the default grid, so the fit
#' # warned (class poemed_grid_boundary). A grid that extends below it and
#' # is spaced more finely moves the choice inside the grid and keeps the
#' # same three mediators.
#' finer <- pe_mediation_linear(d$X, d$Y, d$M,
#'                              lambda_grid = seq(0.01, 10, length.out = 200))
#' attr(finer, "tuning")$lambda_selected
#' attr(finer, "selected_mediators")
#'
#' @export
pe_mediation_linear <- function(X, Y, M, Z = NULL,
                                method = c("Bonferroni", "BH", "BY"),
                                scale = TRUE, error_level = 0.05,
                                report_all_methods = FALSE,
                                drop_constant = FALSE,
                                lambda_grid = seq(0.05, 10, length.out = 100)) {
  method <- .poemed_match_arg(method)
  .validate_error_level(error_level)
  # The logical flags must each be a single TRUE or FALSE; isTRUE() alone
  # would read NA or "yes" as FALSE without a word.
  .validate_flag(scale, "scale")
  .validate_flag(report_all_methods, "report_all_methods")
  # report_all_methods = TRUE screens with all three multiplicity methods (the
  # chosen one stays primary) so their selected sets can be compared.
  methods <- if (isTRUE(report_all_methods))
    unique(c(method, "Bonferroni", "BH", "BY")) else method
  inp <- .pe_validate_inputs(X, Y, M, Z, drop_constant = drop_constant)
  X <- inp$X; Y <- inp$Y; M <- inp$M; Z <- inp$Z
  n <- inp$n; q <- inp$q; p <- ncol(M)

  # One tuning grid serves the full model and the reduced model of the
  # benchmark test, as in the reference implementation's default.
  .validate_lambda_grid(lambda_grid, "lambda_grid")

  # A continuous outcome must vary and should not be binary; if it looks like
  # counts, warn but proceed (the user may genuinely want the linear model).
  # The constant case is named first, since a single value is not binary and
  # the logistic model would not fit it either.
  if (length(unique(Y)) == 1L)
    stop("`Y` is constant, so there is no variation in the outcome to ",
         "mediate. Check the outcome, or the subset of rows it was taken from.",
         call. = FALSE)
  if (length(unique(Y)) <= 2)
    stop("`Y` appears to be binary. Use pe_mediation_logistic() instead.",
         call. = FALSE)
  if (all(Y >= 0) && all(abs(Y - round(Y)) < .Machine$double.eps^0.5))
    warning("`Y` looks like a count variable; consider pe_mediation_poisson().",
            call. = FALSE)

  if (scale) {
    sc <- .pe_scale_data(X, Y, M, Z, center_y = TRUE)
    X <- sc$X; Y <- sc$Y; M <- sc$M; Z <- sc$Z
  }
  S <- Z

  # --- Penalized fit + HBIC tuning for the full model. The grid of lambdas is
  # scored by HBIC and the minimizer (the last, if tied, matching the reference
  # implementation) is taken. The penalized design [M, X, S] does not change
  # across lambdas, so build it once and reuse it for every HBIC evaluation. ---
  MV_full <- if (length(S) == 0) cbind(M, X) else cbind(M, X, S)
  results <- lapply(lambda_grid, .HBIC_calc_linear, xx = X, yy = Y, mm = M,
                    S = S, n_imp = 0, MV = MV_full)
  hbic <- vapply(results, function(r) as.numeric(r$BIC), numeric(1))
  id <- utils::tail(which(hbic == min(hbic)), 1)
  res_full <- results[[id]]
  alpha0_hat <- res_full$alpha0; alpha1_hat <- res_full$alpha1
  alpha2_hat <- res_full$alpha2

  # --- Reduced-model fit (mediators on the outcome with only the intercept or
  # the confounders), used by the inference routine for the direct-effect
  # likelihood-ratio piece. Its design [M, intercept|S] is likewise constant. ---
  if (length(S) == 0) {
    intcpt <- matrix(rep(1, n), ncol = 1)
    MV_red <- cbind(M, intcpt)
    results0 <- lapply(lambda_grid, .HBIC_calc_linear, xx = intcpt,
                       yy = Y, mm = M, n_imp = 0, MV = MV_red)
  } else {
    MV_red <- cbind(M, S)
    results0 <- lapply(lambda_grid, .HBIC_calc_linear, xx = S,
                       yy = Y, mm = M, n_imp = 0, MV = MV_red)
  }
  hbic0 <- vapply(results0, function(r) as.numeric(r$BIC), numeric(1))
  id0 <- utils::tail(which(hbic0 == min(hbic0)), 1)
  alpha0_tld <- results0[[id0]]$alpha0
  alpha2_tld <- results0[[id0]]$alpha1
  tuning <- .pe_tuning_record(lambda_grid, lambda_grid[id], lambda_grid[id0])

  A     <- which(alpha0_hat != 0)    # selected mediators (support of alpha_m)
  A_tld <- which(alpha0_tld != 0)

  if (length(A) == 0L)
    return(.pe_relabel_mediators(
      .pe_empty_output(n, p, q, "continuous (linear)", method,
                       error_level, methods = methods,
                       tuning = tuning,
                       exposure_names = inp$exposure_names),
      inp$kept, inp$mediator_names))
  # An HBIC choice at either end of a grid may not be the criterion's minimum;
  # say so (the empty-fit warning above already covers the empty case).
  .pe_grid_boundary_warning(tuning)

  M_A <- as.matrix(M[, A])
  if (length(A_tld) == 0L) {
    M_tld <- NULL; alpha0_tld_input <- rep(0, length(A))
  } else {
    M_tld <- as.matrix(M[, A_tld]); alpha0_tld_input <- alpha0_tld[A_tld]
  }

  # --- Benchmark Wald inference: S_n and the total indirect effect beta. ---
  hdmm <- .inference_linear(X, Y, M_A, S = S, M_tld = M_tld,
                            alpha0_hat = alpha0_hat[A],
                            alpha0_tld = alpha0_tld_input,
                            alpha1_hat = alpha1_hat, alpha2_hat = alpha2_hat,
                            alpha2_tld = alpha2_tld)
  stat_hdmm <- as.numeric(hdmm$Sn)
  pval_hdmm <- .pe_chisq_pvalue(stat_hdmm, df = q)

  # --- Power enhancement screen(s) and tidy output. ---
  .pe_relabel_mediators(
    .pe_finalize(X, Y, M_A, S, A, stat_hdmm, pval_hdmm, hdmm$beta_hat,
                 hdmm$var_beta, p, n, q, "gaussian", methods, error_level,
                 "continuous (linear)", tuning = tuning,
                       exposure_names = inp$exposure_names),
    inp$kept, inp$mediator_names)
}
