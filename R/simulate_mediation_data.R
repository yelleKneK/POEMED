# Canonical mediator-coefficient (alpha_m) presets per outcome and pattern,
# taken directly from the article's Monte Carlo designs (Section 4.1 for the
# continuous outcome, Section 4.2 for the binary and count outcomes).
# "homogeneous" means all active mediation effects share a sign; "contrasting"
# means the active effects have opposite signs. With the default tau the
# continuous and binary contrasting presets cancel, so the total indirect
# effect beta = Gamma_x alpha_m is zero (the case where a
# total-indirect-effect test is powerless); the count contrasting preset,
# (0, 0, 0, 0.8, -0.7), gives 0.4 * 0.8 - 0.5 * 0.7 = -0.03 per unit of c1.
# The presets place the nonzero coefficients in the first few mediators and
# pad the rest with zeros to length p. Not exported.
#' @keywords internal
#' @noRd
.pe_alpha_m_preset <- function(outcome, pattern, p) {
  nz <- switch(paste(outcome, pattern, sep = "/"),
    "continuous/homogeneous" = c(1, 0.8, 0.6, 0.4, 0.2),
    "continuous/contrasting" = c(1, -0.5, 0.4, -0.3),
    "binary/homogeneous"     = c(3, 1.5, 0, 0, 2),
    "binary/contrasting"     = c(3, -1.5),
    "count/homogeneous"      = c(0.9, 0.8, 0, 0, 0.7),
    "count/contrasting"      = c(0, 0, 0, 0.8, -0.7),
    stop("No preset for outcome/pattern combination.", call. = FALSE))
  if (length(nz) > p)
    stop("`p` is too small for this pattern (need at least ", length(nz),
         " mediators).", call. = FALSE)
  c(nz, rep(0, p - length(nz)))
}

#' Simulate High-Dimensional Mediation Data
#'
#' Generates a data set from the mediation models studied in the article,
#' for a continuous, binary, or count outcome, under the homogeneous or
#' contrasting (heterogeneous) mediation patterns. The mediators carry an
#' AR(1) correlation structure, a handful are truly active, and the rest
#' are null, so the data exercise exactly the situation the power-enhanced
#' tests are designed for. This is the generator behind the package's
#' example data sets and behind [pe_power_curve()].
#'
#' @details
#' One exposure column (or \eqn{q}{q} of them) is drawn standard normal. The
#' mediators follow \eqn{M = X \Gamma_x + Z \Gamma_z + \varepsilon_m}{M = X
#' Gamma_x + Z Gamma_z + eps_m} with \eqn{\varepsilon_m \sim N(0,
#' \Sigma)}{eps_m ~ N(0, Sigma)}, \eqn{\Sigma}{Sigma} the AR(1) covariance
#' [mediation_ar1_cov()] with parameter `rho`. The exposure-on-mediator
#' coefficient is \eqn{\Gamma_x = c_1 \tau}{Gamma_x = c1 tau}, where the
#' base pattern \eqn{\tau}{tau} has small increasing loadings on the first
#' five mediators and random noise loadings on the rest, so `c1` scales the
#' overall exposure-to-mediator signal (and `c1 = 0` gives the global null
#' of no mediation). The outcome is generated from its model with mediator
#' coefficients \eqn{\alpha_m}{alpha_m} (the preset for the chosen
#' `pattern`, or a user-supplied vector) and direct effect \eqn{\alpha_x =
#' c_2}{alpha_x = c2}:
#' \itemize{
#'   \item continuous: \eqn{Y = M\alpha_m + X\alpha_x + Z\alpha_z +
#'     \varepsilon_y}{Y = M alpha_m + X alpha_x + Z alpha_z + eps_y},
#'     \eqn{\varepsilon_y \sim N(0, \sigma_y^2)}{eps_y ~ N(0, sigma_y^2)};
#'   \item binary: \eqn{Y \sim \mathrm{Bernoulli}(\mathrm{logit}^{-1}(
#'     M\alpha_m + X\alpha_x + Z\alpha_z))}{Y ~ Bernoulli(plogis(M alpha_m
#'     + X alpha_x + Z alpha_z))};
#'   \item count: \eqn{Y \sim \mathrm{Poisson}(\exp(\eta))}{Y ~
#'     Poisson(exp(eta))}, with each value of the log-mean
#'     \eqn{\eta}{eta} clamped to the interval from -5 to 5 (a value outside
#'     it is set to the nearer endpoint) to guard against overflow.
#' }
#'
#' Every design argument is checked before any random number is drawn: a
#' coefficient vector or matrix of the wrong length or shape, or one holding
#' a missing or infinite value, stops with a message naming it.
#'
#' @param n Number of observations.
#' @param p Number of candidate mediators.
#' @param outcome Outcome type: `"continuous"`, `"binary"`, or `"count"`.
#' @param pattern Mediation pattern: `"homogeneous"` (active effects share
#'   a sign) or `"contrasting"` (active effects of opposite sign).
#'   `"heterogeneous"` is accepted as a synonym for `"contrasting"`. Ignored
#'   when `alpha_m` is supplied directly. The presets are the article's. With
#'   the default `tau`, the contrasting presets for the continuous and binary
#'   outcomes cancel to a total indirect effect of zero, while the count
#'   preset, (0, 0, 0, 0.8, -0.7), gives a total indirect effect of -0.03
#'   times `c1`.
#' @param c1 Scale on the exposure-on-mediator coefficients \eqn{\Gamma_x
#'   = c_1 \tau}{Gamma_x = c1 tau}. `c1 = 0` makes every mediator inactive
#'   (the global null). Default 1.
#' @param c2 The direct effect \eqn{\alpha_x}{alpha_x} of each exposure on
#'   the outcome. Default 0.5.
#' @param q Number of exposures. Default 1. For `q > 1` you must supply
#'   `tau` as a `q` by `p` matrix.
#' @param d Number of confounders. Default 0 (none). When positive,
#'   confounders are drawn standard normal and enter both models through
#'   `alpha_z` and `Gamma_z`.
#' @param alpha_m Optional numeric vector of `p` finite mediator-on-outcome
#'   coefficients, overriding the `pattern` preset.
#' @param tau Optional base pattern for \eqn{\Gamma_x}{Gamma_x}, finite
#'   numbers. By default it is the article's (0.1, 0.2, 0.3, 0.4, 0.5,
#'   noise), with the noise loadings drawn from N(0, 0.5^2), for `q = 1`.
#'   Required when `q > 1`, as a `q` by `p` matrix (one row per exposure)
#'   or a vector of its values in column order; when `q = 1` a vector of
#'   length `p` serves.
#' @param rho AR(1) correlation parameter for the mediator noise. Default
#'   0.5.
#' @param sigma_y Outcome noise standard deviation for a continuous
#'   outcome. Default 0.5.
#' @param alpha_z,Gamma_z Optional confounder coefficients: a vector of `d`
#'   finite numbers and a `d` by `p` matrix (when `d = 1`, a vector of
#'   length `p` serves for `Gamma_z`). Each defaults to zeros when `d > 0`.
#'   They are used only when `d > 0`; supplying either with `d = 0`, where
#'   it would be ignored, is an error.
#' @param seed Optional integer seed. When supplied it is set for the
#'   duration of the call and the caller's random number generator state is
#'   restored on exit. `NULL` (the default) sets no seed: the draws come from
#'   the session's random number stream, which the call advances as any
#'   random function does.
#'
#' @return A list with the simulated data and the truth used to generate
#'   it: `X` (`n` by `q`), `M` (`n` by `p`), `Y` (length `n`), `Z` (`n` by
#'   `d`, or `NULL`), the coefficient values `alpha_m`, `Gamma_x`,
#'   `alpha_x`, the total indirect effect `beta` (`Gamma_x` times
#'   `alpha_m`), and `active_mediators` (the indices of the truly active
#'   mediators, those nonzero in both paths), together with `n`, `p`, `q`,
#'   `outcome`, and `pattern`.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @references
#' Yu, X., & Kelley, K. (in press). Power Enhancement in
#' High-Dimensional Heterogeneous Mediation Analysis. \emph{Journal of the
#' American Statistical Association}.
#'
#' @seealso [pe_mediation()] to analyze the data, [pe_power_curve()] to
#'   sweep `c1`, [mediation_ar1_cov()] for the mediator covariance.
#'
#' @family mediation simulation
#'
#' @examples
#' set.seed(113)
#' d <- simulate_mediation_data(n = 100, p = 50, outcome = "continuous",
#'                              pattern = "contrasting")
#' dim(d$M)
#' d$active_mediators       # truly active mediators
#' round(d$beta, 10)        # total indirect effect, zero up to rounding error
#'
#' @export
simulate_mediation_data <- function(n, p,
                                     outcome = c("continuous", "binary", "count"),
                                     pattern = c("homogeneous", "contrasting",
                                                 "heterogeneous"),
                                     c1 = 1, c2 = 0.5, q = 1L, d = 0L,
                                     alpha_m = NULL, tau = NULL, rho = 0.5,
                                     sigma_y = 0.5, alpha_z = NULL,
                                     Gamma_z = NULL, seed = NULL) {
  # Validate the whole design, and resolve the presets, before any draw.
  des <- .pe_sim_design(n = n, p = p, outcome = outcome, pattern = pattern,
                        c1 = c1, c2 = c2, q = q, d = d, alpha_m = alpha_m,
                        tau = tau, rho = rho, sigma_y = sigma_y,
                        alpha_z = alpha_z, Gamma_z = Gamma_z)
  outcome <- des$outcome
  pattern <- des$pattern
  alpha_m <- des$alpha_m

  # Local seeding that restores the caller's RNG state on exit (see poemed_seed.R).
  .poemed_local_seed(seed)

  # Base exposure-on-mediator pattern tau (q-by-p), scaled by c1 to form
  # Gamma_x. The default single-exposure tau loads 0.1..0.5 on the first five
  # mediators and N(0, 0.5^2) noise on the rest, following Yu and Kelley (in press).
  if (is.null(tau)) {
    tau <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5,
                    stats::rnorm(p - 5, 0, 0.5)), nrow = 1L)
  } else {
    tau <- matrix(tau, nrow = q)
  }
  Gamma_x <- c1 * tau                          # q-by-p
  alpha_x <- rep(c2, q)                         # direct effect per exposure

  # Exposures, correlated mediator noise, and mediators.
  X <- matrix(stats::rnorm(n * q), nrow = n, ncol = q)
  Sigma_m <- mediation_ar1_cov(p, rho = rho)
  eps_m <- matrix(stats::rnorm(n * p), nrow = n, ncol = p) %*% chol(Sigma_m)
  M <- X %*% Gamma_x + eps_m

  # Optional confounders enter both the mediator and outcome models.
  Z <- NULL
  if (d > 0L) {
    Z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
    Gamma_z <- des$Gamma_z
    alpha_z <- des$alpha_z
    M <- M + Z %*% Gamma_z
  }

  # Linear predictor for the outcome model.
  eta <- as.numeric(M %*% alpha_m + X %*% alpha_x)
  if (d > 0L) eta <- eta + as.numeric(Z %*% alpha_z)

  Y <- switch(outcome,
    continuous = eta + stats::rnorm(n, 0, sigma_y),
    binary     = stats::rbinom(n, 1, stats::plogis(eta)),
    # Clamp each value of the log-mean to [-5, 5] before exponentiating, so a
    # few large linear predictors do not produce overflowing counts.
    count      = stats::rpois(n, exp(pmax(pmin(eta, 5), -5))))

  beta <- as.numeric(Gamma_x %*% alpha_m)       # total indirect effect
  # Truly active mediators: nonzero on BOTH paths (outcome and exposure).
  active <- which(alpha_m != 0 & apply(Gamma_x != 0, 2L, any))

  list(X = X, M = M, Y = Y, Z = Z, alpha_m = alpha_m, Gamma_x = Gamma_x,
       alpha_x = alpha_x, beta = beta, active_mediators = active,
       n = n, p = p, q = q, outcome = outcome, pattern = pattern)
}

# The design half of simulate_mediation_data(): validate every argument by
# name and resolve the pieces that need no random numbers (the outcome and
# pattern, the alpha_m preset, the zero confounder coefficients). It draws
# nothing, so the study functions also call it once before their replications
# to stop a bad design with the simulator's own message (inside a forked
# worker the message would be lost). Its arguments are those of
# simulate_mediation_data() without `seed`, with the same defaults; a test
# keeps the two in step. Returns a list with `outcome`, `pattern`, `alpha_m`,
# `alpha_z`, and `Gamma_z`. Not exported.
#' @keywords internal
#' @noRd
.pe_sim_design <- function(n, p,
                           outcome = c("continuous", "binary", "count"),
                           pattern = c("homogeneous", "contrasting",
                                       "heterogeneous"),
                           c1 = 1, c2 = 0.5, q = 1L, d = 0L,
                           alpha_m = NULL, tau = NULL, rho = 0.5,
                           sigma_y = 0.5, alpha_z = NULL, Gamma_z = NULL) {
  outcome <- .pe_match_choice(outcome, c("continuous", "binary", "count"),
                              "outcome")
  pattern <- .pe_match_choice(pattern, c("homogeneous", "contrasting",
                                         "heterogeneous"), "pattern")
  if (pattern == "heterogeneous") pattern <- "contrasting"
  .validate_sim_design(n, p, q, d, c1, c2, rho, sigma_y)

  finite_numbers <- function(x) is.numeric(x) && all(is.finite(x))

  # Mediator-on-outcome coefficients: a user vector or the pattern preset.
  if (is.null(alpha_m)) {
    alpha_m <- .pe_alpha_m_preset(outcome, pattern, p)
  } else {
    if (!finite_numbers(alpha_m))
      stop("`alpha_m` must be a numeric vector of finite values, one ",
           "coefficient per mediator.", call. = FALSE)
    if (length(alpha_m) != p)
      stop(sprintf("`alpha_m` must have length `p` (%d), one coefficient per mediator; it has length %d.",
                   as.integer(p), length(alpha_m)), call. = FALSE)
  }

  # The exposure-on-mediator base pattern: the default needs q = 1 and at
  # least five mediators. A supplied tau is a q-by-p matrix, or a vector of
  # q * p values that matrix(tau, nrow = q) reads into it by column; with a
  # single exposure any shape holding p values serves. A matrix of another
  # shape (a p-by-q one, say) would be silently rearranged, so it stops.
  if (is.null(tau)) {
    if (q != 1L)
      stop("For `q > 1` you must supply `tau` as a q-by-p matrix.",
           call. = FALSE)
    if (p < 5L) stop("`p` must be at least 5 for the default `tau`.",
                     call. = FALSE)
  } else {
    shape_ok <- if (q == 1L) {
      length(tau) == p
    } else if (is.matrix(tau)) {
      identical(as.integer(dim(tau)), as.integer(c(q, p)))
    } else {
      length(tau) == q * p
    }
    if (!finite_numbers(tau) || !shape_ok)
      stop(sprintf(paste0("`tau` must be a %d-by-%d matrix of finite numbers ",
                          "(one row per exposure, one column per mediator), ",
                          "or a vector of its %d values in column order."),
                   as.integer(q), as.integer(p), as.integer(q * p)),
           call. = FALSE)
  }

  # Confounder coefficients: meaningful only with confounders (d > 0), where
  # each defaults to zeros; supplied with d = 0 they would be silently
  # ignored, so that stops instead.
  if (d == 0L) {
    given <- c(if (length(alpha_z)) "`alpha_z`", if (length(Gamma_z)) "`Gamma_z`")
    if (length(given))
      stop(paste(given, collapse = " and "),
           if (length(given) == 1L) " describes" else " describe",
           " confounders, but `d = 0` asks for none, so ",
           if (length(given) == 1L) "it" else "they",
           " would be ignored. Set `d` to the number of confounders, or ",
           "leave ", paste(given, collapse = " and "), " at NULL.",
           call. = FALSE)
  } else {
    if (is.null(alpha_z)) {
      alpha_z <- rep(0, d)
    } else if (!finite_numbers(alpha_z) || length(alpha_z) != d) {
      stop(sprintf("`alpha_z` must be a numeric vector of %d finite values, one per confounder (`d` = %d).",
                   as.integer(d), as.integer(d)), call. = FALSE)
    }
    if (is.null(Gamma_z)) {
      Gamma_z <- matrix(0, nrow = d, ncol = p)
    } else {
      shape_ok <- if (is.matrix(Gamma_z))
        identical(as.integer(dim(Gamma_z)), as.integer(c(d, p)))
      else d == 1L && length(Gamma_z) == p
      if (!finite_numbers(Gamma_z) || !shape_ok)
        stop(sprintf(paste0("`Gamma_z` must be a %d-by-%d matrix of finite ",
                            "numbers (one row per confounder, one column per ",
                            "mediator)%s."),
                     as.integer(d), as.integer(p),
                     if (d == 1L) ", or a numeric vector of length `p`" else ""),
             call. = FALSE)
    }
  }

  list(outcome = outcome, pattern = pattern, alpha_m = alpha_m,
       alpha_z = alpha_z, Gamma_z = Gamma_z)
}
