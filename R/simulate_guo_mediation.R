#' Simulate the Real-Data-Motivated Heterogeneous Mediation Setting
#'
#' Generates a linear mediation data set from the article's
#' \emph{real-data-motivated} heterogeneous setting, whose coefficients are
#' calibrated to the DNA-methylation case study of Guo et al. (2022). Unlike
#' the homogeneous and contrasting presets of [simulate_mediation_data()],
#' this setting has eleven active mediators with a genuine \strong{mix of
#' positive and negative} effects that neither all agree in sign (homogeneous)
#' nor exactly cancel (contrasting): the messy middle ground that real data
#' tends to present. The calibration constants are the [guo_calibration]
#' object; the total indirect effect is `guo_calibration$beta_per_c1` times
#' \eqn{c_1}{c1}, about \eqn{-1.597 c_1}{-1.597 c1}. (The article reports
#' \eqn{-1.5977 c_1}{-1.5977 c1}, from the full-precision Guo coefficients;
#' the values shipped here are those coefficients rounded to three decimals.)
#'
#' The mediator dimension (\eqn{p = 1008}{p = 1008}), the eleven active
#' loci, and the coefficient vectors come from Guo et al. (2022) and are
#' fixed. The exposure and confounders, which the article took from the
#' actual case-study data, are simulated here (standard normal) so the
#' function is self-contained. The mediator and outcome coefficients, the
#' part of the calibration that matters scientifically, are the values
#' printed in the supplement of Yu and Kelley (in press), where the Guo
#' estimates appear to three decimals.
#'
#' @param c1 Signal-strength scale for the exposure-mediator coefficients
#'   (`Gamma_x = c1 * Gamma_x_G`), a single finite number. At `c1 = 0` there
#'   is no mediation (the Type I error setting); the article sweeps `c1`
#'   over `c(0, +/- 0.1, ..., +/- 1)`.
#' @param c2 Direct effect (exposure-outcome coefficient), a single finite
#'   number. Default 0.5.
#' @param confounders `TRUE` or `FALSE`; if `TRUE`, include the calibrated
#'   confounders in the data-generating process (the article's "with
#'   confounders" scenario). Default `FALSE`.
#' @param alpha_m Outcome-mediator coefficients at the eleven active loci:
#'   either a numeric vector of 11 finite values, or one of the strings
#'   `"setting0"` (the default, the values 1, 0.9, 0.8, -0.9, -0.8, -0.7,
#'   0.6, 0.5, 0.4, 0.3, and 0.2), `"homogeneous_like"`, or
#'   `"contrasting_like"` (the two alternative sets the article also
#'   studies, which remain heterogeneous because they are paired with the
#'   scattered `Gamma_x`).
#' @param n Number of observations, a single whole number of at least 2.
#'   Default `NULL`, the case-study size, 85.
#' @param seed Optional integer seed. When supplied it is set for the
#'   duration of the call and the caller's random number generator state is
#'   restored on exit. `NULL` (the default) sets no seed: the draws come from
#'   the session's random number stream, which the call advances as any
#'   random function does.
#'
#' @return A list with the same shape as [simulate_mediation_data()]: the
#'   exposure `X` (`n` by 1), mediators `M` (`n` by 1008), outcome `Y`,
#'   confounders `Z` (`n` by 8, or `NULL`), the full coefficient vectors
#'   `alpha_m` and `Gamma_x`, the total indirect effect `beta`, the
#'   `active_mediators` (the loci nonzero in both paths: all eleven under the
#'   default `alpha_m`, the first five under `"homogeneous_like"`, the first
#'   four under `"contrasting_like"`, and none when `c1 = 0`), and the
#'   dimensions.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @references
#' Yu, X., & Kelley, K. (in press). Power Enhancement in High-Dimensional
#' Heterogeneous Mediation Analysis. \emph{Journal of the American Statistical
#' Association}.
#'
#' Guo, X., Li, R., Liu, J., & Zeng, M. (2022). High-dimensional mediation
#' analysis for selecting DNA methylation loci mediating childhood trauma and
#' cortisol stress reactivity. \emph{Journal of the American Statistical
#' Association, 117}(539), 1110--1121.
#' \doi{10.1080/01621459.2022.2053136}
#'
#' @seealso [simulate_mediation_data()], [guo_calibration],
#'   [pe_power_curve()].
#'
#' @family mediation simulation
#'
#' @examples
#' set.seed(113)
#' d <- simulate_guo_mediation(c1 = 0.5)
#' dim(d$M)               # 85 x 1008
#' d$active_mediators     # the 11 active loci
#' round(d$beta, 4)       # beta_per_c1 * 0.5, i.e. -0.7985
#' # At the case study's n = 85 and c1 = 0.5 the signal is faint for either
#' # test. With c1 = 1 and n = 150 both tests reject the global null, and the
#' # power-enhanced test also selects seven of the eleven active loci. The
#' # fit searches a 20-value grid to keep the example fast; the default grid
#' # selects the same seven loci and takes five times as long.
#' d2 <- simulate_guo_mediation(c1 = 1, n = 150, seed = 113)
#' fit <- pe_mediation(d2$X, d2$Y, d2$M, outcome = "continuous",
#'                     lambda_grid = seq(0.05, 2, length.out = 20))
#' fit[fit$term %in% c("pval_hdmm", "pval_pe", "n_selected_mediators"), ]
#'
#' @export
simulate_guo_mediation <- function(c1, c2 = 0.5, confounders = FALSE,
                                   alpha_m = "setting0", n = NULL,
                                   seed = NULL) {
  cal <- guo_calibration
  p <- cal$p; loc <- cal$locations
  if (is.null(n)) n <- cal$n
  # The same checks, and messages, as simulate_mediation_data() for the
  # arguments the two simulators share.
  if (!is.numeric(c1) || length(c1) != 1L || !is.finite(c1))
    stop("`c1` must be a single finite number.", call. = FALSE)
  if (!is.numeric(c2) || length(c2) != 1L || !is.finite(c2))
    stop("`c2` must be a single finite number.", call. = FALSE)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 2 ||
      n != round(n))
    stop("`n` must be a single whole number of at least 2.", call. = FALSE)
  .pe_check_flag(confounders, "confounders")

  # Resolve the 11 active-locus outcome-mediator coefficients.
  am11 <- if (is.character(alpha_m)) {
    switch(.pe_match_choice(alpha_m, c("setting0", "homogeneous_like",
                                       "contrasting_like"), "alpha_m"),
           setting0 = cal$alpha_m,
           homogeneous_like = cal$alpha_m_variants$homogeneous_like[1:11],
           contrasting_like = cal$alpha_m_variants$contrasting_like[1:11])
  } else {
    if (!is.numeric(alpha_m) || !all(is.finite(alpha_m)))
      stop("A numeric `alpha_m` must hold finite values, one per active ",
           "locus.", call. = FALSE)
    if (length(alpha_m) != length(loc))
      stop("Numeric `alpha_m` must have length ", length(loc), ".",
           call. = FALSE)
    alpha_m
  }

  .poemed_local_seed(seed)

  # Full-length coefficient vectors (zero off the active loci).
  alpha_m_full <- numeric(p); alpha_m_full[loc] <- am11
  Gamma_x_full <- numeric(p); Gamma_x_full[loc] <- cal$Gamma_x
  Gamma_x_mat <- matrix(c1 * Gamma_x_full, nrow = 1L)        # 1 x p

  X <- matrix(stats::rnorm(n), nrow = n, ncol = 1L)
  Sigma_m <- mediation_ar1_cov(p, rho = 0.5)
  eps_m <- matrix(stats::rnorm(n * p), nrow = n, ncol = p) %*% chol(Sigma_m)
  M <- X %*% Gamma_x_mat + eps_m

  Z_fit <- NULL; eta_z <- 0
  if (confounders) {
    # Nine confounders: an intercept column plus eight covariates. The
    # intercept's constant contribution is absorbed by centering when fitting,
    # so the returned Z holds the eight varying covariates.
    Zd <- cbind(1, matrix(stats::rnorm(n * 8L), nrow = n, ncol = 8L))
    M[, loc] <- M[, loc] + Zd %*% t(cal$Gamma_z)
    eta_z <- as.numeric(Zd %*% cal$alpha_z)
    Z_fit <- Zd[, -1L, drop = FALSE]
  }

  eta <- as.numeric(M %*% alpha_m_full + X * c2 + eta_z)
  Y <- eta + stats::rnorm(n, 0, 0.5)

  beta <- c1 * sum(cal$Gamma_x * am11)
  # Truly active loci: nonzero on BOTH paths, as in simulate_mediation_data().
  # Every calibrated Gamma_x is nonzero, so with c1 != 0 these are the loci
  # whose outcome coefficient is nonzero.
  active <- if (c1 != 0) loc[am11 != 0 & cal$Gamma_x != 0] else integer(0)

  list(X = X, M = M, Y = Y, Z = Z_fit, alpha_m = alpha_m_full,
       Gamma_x = Gamma_x_mat, beta = beta, active_mediators = active,
       n = n, p = p, q = 1L, outcome = "continuous", pattern = "guo")
}
