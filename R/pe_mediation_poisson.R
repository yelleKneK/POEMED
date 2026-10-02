#' Power-Enhanced Mediation Test for a Count Outcome
#'
#' Tests the global null hypothesis of no active mediator for a count
#' outcome modeled with Poisson (log-link) regression, reporting both
#' the benchmark Wald test on the total indirect effect and the
#' power-enhanced (PE) test of Yu and Kelley (in press). This is the
#' count-outcome worker behind [pe_mediation()].
#'
#' @details
#' The outcome follows a Poisson mediation model,
#' \eqn{Y \mid M, X, Z \sim \mathrm{Poisson}(\exp(\alpha_m' M +
#' \alpha_x' X + \alpha_z' Z))}{Y | M, X, Z ~ Poisson(exp(alpha_m' M +
#' alpha_x' X + alpha_z' Z))},
#' with a linear mediator model as in [pe_mediation_linear()]. Estimation,
#' tuning, and the benchmark Wald test follow the generalized-outcome
#' framework of Guo et al. (2024): the mediator coefficients are estimated by
#' partial penalized likelihood with a SCAD penalty (Fan and Li, 2001),
#' leaving the exposures and confounders unpenalized, with the tuning
#' parameter chosen by the high-dimensional BIC of Wang, Kim, and Li (2013).
#' The penalized Poisson fit is a glmnet lasso fit followed by one weighted
#' lasso step on a quadratic approximation of the Poisson likelihood at the
#' lasso estimate, with the SCAD derivatives at the lasso coefficients as
#' weights: the one-step local linear approximation of Zou and Li (2008),
#' started from the lasso as in Fan, Xue, and Zou (2014). The
#' power-enhanced test of Yu and Kelley (in press), built on the power
#' enhancement principle of Fan, Liao, and Yao (2015), adds the component
#' \eqn{J_m}{J_m} (see [pe_mediation()]). The total indirect effect
#' \eqn{\beta = \Gamma_x \alpha_m}{beta = Gamma_x alpha_m} is on the log-mean
#' (link) scale, not the count scale: with `scale = TRUE` it is the change in
#' the log mean count per standard deviation of the exposure.
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
#' @seealso [pe_mediation()], [pe_mediation_linear()],
#'   [pe_mediation_logistic()].
#'
#' @family mediation tests
#'
#' @examples
#' set.seed(113)
#' d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
#'                              outcome = "count", c1 = 0.5, c2 = 0.4)
#' pe_mediation_poisson(d$X, d$Y, d$M)
#'
#' @export
pe_mediation_poisson <- function(X, Y, M, Z = NULL,
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
  phi0 <- 1                          # known dispersion for the Poisson family
  .validate_lambda_grid(lambda_grid, "lambda_grid")

  # The outcome must be a nonnegative integer count with some variation.
  if (any(Y < 0, na.rm = TRUE))
    stop("`Y` must be nonnegative for the Poisson model.", call. = FALSE)
  if (any(abs(Y - round(Y)) > sqrt(.Machine$double.eps), na.rm = TRUE))
    stop("`Y` must contain integer counts for the Poisson model.", call. = FALSE)
  if (length(unique(Y[!is.na(Y)])) == 1L)
    stop("`Y` has no variation; the Poisson model is not identifiable.",
         call. = FALSE)

  if (scale) {
    sc <- .pe_scale_data(X, Y, M, Z, center_y = FALSE)
    X <- sc$X; Y <- sc$Y; M <- sc$M; Z <- sc$Z
  }
  S <- Z

  # --- Penalized Poisson fit + HBIC tuning. (The Poisson tie-break takes the
  # first minimizer, matching the reference implementation.) ---
  ngrid <- length(lambda_grid)
  hbic <- numeric(ngrid)
  alpha_rcd <- matrix(0, nrow = ngrid, ncol = p + q + s)
  for (j in seq_len(ngrid)) {
    fit <- .HBIC_poisson(X, Y, M, lambda_grid[j], S)
    hbic[j] <- fit$HBIC
    alpha_rcd[j, ] <- fit$alpha_hat
  }
  id <- which(hbic == min(hbic))[1]
  alpha_hat <- alpha_rcd[id, ]
  alpha0_hat <- alpha_hat[1:p]
  alpha1_hat <- alpha_hat[(p + 1):(p + q)]
  alpha2_hat <- if (s > 0) alpha_hat[(p + q + 1):(p + q + s)] else NULL
  A <- which(alpha0_hat != 0)
  tuning <- .pe_tuning_record(lambda_grid, lambda_grid[id])

  if (length(A) == 0L)
    return(.pe_relabel_mediators(
      .pe_empty_output(n, p, q, "count (Poisson)", method,
                       error_level, methods = methods,
                       tuning = tuning,
                       exposure_names = inp$exposure_names),
      inp$kept, inp$mediator_names))
  .pe_grid_boundary_warning(tuning)

  M_A <- as.matrix(M[, A])

  # --- Benchmark Wald inference for the Poisson model. ---
  hdmm <- .Testing_poisson(X, Y, M_A, S, phi0, alpha0_hat = alpha0_hat[A],
                           alpha1_hat = alpha1_hat, alpha2_hat = alpha2_hat)
  stat_hdmm <- as.numeric(hdmm$Sn)
  pval_hdmm <- .pe_chisq_pvalue(stat_hdmm, df = q)

  # --- Power enhancement component (outcome glm is Poisson). ---
  # --- Power enhancement screen(s) and tidy output. ---
  .pe_relabel_mediators(
    .pe_finalize(X, Y, M_A, S, A, stat_hdmm, pval_hdmm, hdmm$beta_hat,
                 hdmm$var_beta, p, n, q, "poisson", methods, error_level,
                 "count (Poisson)", tuning = tuning,
                       exposure_names = inp$exposure_names),
    inp$kept, inp$mediator_names)
}
