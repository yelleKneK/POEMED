#' Power-Enhanced Mediation Test for a Binary Outcome
#'
#' Tests the global null hypothesis of no active mediator for a binary
#' (0/1) outcome modeled with logistic regression, reporting both the
#' benchmark Wald test on the total indirect effect and the
#' power-enhanced (PE) test of Yu and Kelley (in press). This is the
#' binary-outcome worker behind [pe_mediation()].
#'
#' @details
#' The outcome follows a logistic mediation model,
#' \eqn{\mathrm{logit}\,P(Y = 1 \mid M, X, Z) = \alpha_m' M + \alpha_x' X
#' + \alpha_z' Z}{logit P(Y = 1 | M, X, Z) = alpha_m' M + alpha_x' X +
#' alpha_z' Z},
#' with a linear mediator model as in [pe_mediation_linear()]. Following Guo
#' et al. (2024), the mediator coefficients are estimated by partial
#' penalized likelihood, leaving the exposures and confounders unpenalized,
#' with the tuning parameter chosen by the high-dimensional BIC of Wang, Kim,
#' and Li (2013). The benchmark Wald test on the total indirect effect
#' \eqn{\beta = \Gamma_x \alpha_m}{beta = Gamma_x alpha_m} is also that of
#' Guo et al. (2024). The power-enhanced test of Yu and Kelley (in press),
#' built on the power enhancement principle of Fan, Liao, and Yao (2015),
#' adds the component \eqn{J_m}{J_m}, formed on the logit (link) scale
#' exactly as in the continuous case (see [pe_mediation()]). The total
#' indirect effect here is on the log-odds scale, not the probability scale:
#' with `scale = TRUE` it is the change in the log odds per standard
#' deviation of the exposure.
#'
#' The penalty is SCAD (Fan and Li, 2001) with an adaptive rescaling, not
#' the plain SCAD penalty. The fit is an iteratively reweighted coordinate
#' descent. Each coordinate update applies the SCAD thresholds to the
#' mediator's weighted working score and only then divides by its
#' curvature \eqn{v_j}{v_j} (the mean of the squared standardized mediator
#' weighted by the working weights, at most 1/4 for a logistic fit). The
#' fit is therefore a stationary point of the penalized likelihood with
#' penalty
#' \eqn{\sum_j p_\lambda(v_j |\alpha_{m,j}|) / v_j}{
#' sum_j p_lambda(v_j |alpha_m[j]|) / v_j},
#' not with the plain SCAD penalty
#' \eqn{\sum_j p_\lambda(|\alpha_{m,j}|)}{sum_j p_lambda(|alpha_m[j]|)}. The
#' condition for a coefficient to be zero is the same as under plain SCAD
#' (its score is at most \eqn{\lambda}{lambda} in absolute value), but a
#' nonzero coefficient escapes shrinkage only beyond
#' \eqn{\gamma \lambda / v_j}{gamma lambda / v_j} rather than
#' \eqn{\gamma \lambda}{gamma lambda} (\eqn{\gamma = 3.7}{gamma = 3.7}), so
#' coefficients of moderate size are shrunk more than plain SCAD would shrink
#' them.
#'
#' @inheritParams pe_mediation
#'
#' @return A tidy `data.frame` of class `poemed_tbl` with the same rows
#'   and attributes as [pe_mediation_linear()], except that the
#'   `"tuning"` list has no reduced-model entry.
#'
#' @references
#' Yu, X., & Kelley, K. (in press). Power Enhancement in
#' High-Dimensional Heterogeneous Mediation Analysis. \emph{Journal of the
#' American Statistical Association}.
#'
#' Guo, X., Li, R., Liu, J., & Zeng, M. (2024). Estimations and tests
#' for generalized mediation models with high-dimensional potential
#' mediators. \emph{Journal of Business & Economic Statistics, 42}(1),
#' 243--256. \doi{10.1080/07350015.2023.2174548}
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
#' @seealso [pe_mediation()], [pe_mediation_linear()],
#'   [pe_mediation_poisson()].
#'
#' @family mediation tests
#'
#' @examples
#' set.seed(113)
#' # Contrasting mediation with a binary outcome: the benchmark Wald test
#' # does not reject, the PE test does and selects one of the two active
#' # mediators.
#' d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
#'                              outcome = "binary", c1 = 1, c2 = 1)
#' pe_mediation_logistic(d$X, d$Y, d$M)
#' # HBIC chose 0.05, the smallest value of the default grid, so the fit
#' # warned (class poemed_grid_boundary); a grid that extends below it and
#' # is spaced more finely, seq(0.01, 10, length.out = 200), selects the
#' # same mediator.
#'
#' @export
pe_mediation_logistic <- function(X, Y, M, Z = NULL,
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
  # report_all_methods = TRUE screens with all three multiplicity methods
  # (the chosen one stays primary) so their selected sets can be compared.
  methods <- if (isTRUE(report_all_methods))
    unique(c(method, "Bonferroni", "BH", "BY")) else method
  inp <- .pe_validate_inputs(X, Y, M, Z, drop_constant = drop_constant)
  X <- inp$X; Y <- inp$Y; M <- inp$M; Z <- inp$Z
  n <- inp$n; q <- inp$q; s <- inp$s; p <- ncol(M)
  phi0 <- 1                          # known dispersion for the binomial family
  .validate_lambda_grid(lambda_grid, "lambda_grid")

  # The outcome must be binary and coded 0/1 for the logistic model.
  y_unique <- unique(Y[!is.na(Y)])
  if (length(y_unique) == 1L)
    stop("`Y` is constant (every value is ", y_unique, "); the logistic model ",
         "needs both 0s and 1s in the outcome.", call. = FALSE)
  if (length(y_unique) != 2L)
    stop("`Y` must be binary (two distinct values) for the logistic model.",
         call. = FALSE)
  if (!all(y_unique %in% c(0, 1)))
    stop("`Y` must be coded as 0/1 for the logistic model. Detected values: ",
         paste(y_unique, collapse = ", "), ".", call. = FALSE)

  if (scale) {
    sc <- .pe_scale_data(X, Y, M, Z, center_y = FALSE)
    X <- sc$X; Y <- sc$Y; M <- sc$M; Z <- sc$Z
  }
  S <- Z

  # --- Penalized partial-likelihood fit + HBIC tuning. The penalty factor w
  # penalizes only the mediator block (first p columns); exposures and
  # confounders are always retained. ---
  w <- rep(0, p + q + s); w[1:p] <- 1
  ngrid <- length(lambda_grid)
  hbic <- numeric(ngrid)
  alpha_rcd <- matrix(0, nrow = ngrid, ncol = p + q + s)
  for (j in seq_len(ngrid)) {
    fit <- .HBIC_bino(X = X, Y = Y, M = M, S = S, w = w, lamb = lambda_grid[j])
    hbic[j] <- fit$HBIC
    alpha_rcd[j, ] <- fit$alpha
  }
  id <- utils::tail(which(hbic == min(hbic)), 1)
  alpha_hat <- alpha_rcd[id, ]
  alpha0_hat <- alpha_hat[1:p]
  alpha1_hat <- alpha_hat[(p + 1):(p + q)]
  alpha2_hat <- if (s > 0) alpha_hat[(p + q + 1):(p + q + s)] else NULL
  A <- which(alpha0_hat != 0)
  tuning <- .pe_tuning_record(lambda_grid, lambda_grid[id])

  if (length(A) == 0L)
    return(.pe_relabel_mediators(
      .pe_empty_output(n, p, q, "binary (logistic)", method,
                       error_level, methods = methods,
                       tuning = tuning,
                       exposure_names = inp$exposure_names),
      inp$kept, inp$mediator_names))
  .pe_grid_boundary_warning(tuning)

  M_A <- as.matrix(M[, A])

  # --- Benchmark Wald inference for the logistic model. ---
  hdmm <- .Testing_bino(X, Y, M_A, S, phi0, a0_hat = alpha0_hat[A],
                        a1_hat = alpha1_hat, a2_hat = alpha2_hat)
  stat_hdmm <- as.numeric(hdmm$Sn)
  pval_hdmm <- .pe_chisq_pvalue(stat_hdmm, df = q)

  # --- Power enhancement component (outcome glm is binomial). ---
  # --- Power enhancement screen(s) and tidy output. ---
  .pe_relabel_mediators(
    .pe_finalize(X, Y, M_A, S, A, stat_hdmm, pval_hdmm, hdmm$beta_hat,
                 hdmm$var_beta, p, n, q, "binomial", methods, error_level,
                 "binary (logistic)", tuning = tuning,
                       exposure_names = inp$exposure_names),
    inp$kept, inp$mediator_names)
}
