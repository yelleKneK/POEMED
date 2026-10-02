# Internal utilities shared by the POEMED testing and simulation functions.
# None are exported. They concentrate the argument coercion, validation, and
# preprocessing that every pe_mediation_*() worker needs, so the workers
# differ only in the parts that are genuinely family-specific (the penalized
# fit and the Wald inference). The reference PEmediation implementation
# repeated this logic in each worker; centralizing it here removes that
# duplication and gives one place to test the input contract.

# Coerce X / M / Z to a numeric double matrix. A bare numeric vector becomes a
# one-column matrix (the q = 1 single-exposure case, which is the common one),
# a data.frame is matricized, and storage mode is forced to double so the
# downstream linear algebra never silently promotes from integer. NULL passes
# through as NULL so an absent confounder matrix Z stays absent.
#' @keywords internal
#' @noRd
.as_numeric_matrix <- function(x, name) {
  if (is.null(x)) return(NULL)
  if (is.vector(x) && !is.list(x)) x <- matrix(x, ncol = 1L)
  if (!is.matrix(x) && !is.data.frame(x))
    stop(sprintf("`%s` must be a numeric matrix, data.frame, or numeric vector.",
                 name), call. = FALSE)
  x <- as.matrix(x)
  # A zero-column matrix (from an empty data.frame, which as.matrix() makes
  # logical) skips the type check, so the caller's column-count guard can say
  # what is missing instead of calling it non-numeric.
  if (ncol(x) > 0L && !is.numeric(x))
    stop(sprintf("`%s` must be numeric.", name), call. = FALSE)
  storage.mode(x) <- "double"
  x
}

# Coerce the outcome Y to a numeric vector. A one-column matrix is flattened
# (a convenience, since Y is naturally a column), anything list-like is
# rejected, and the result is cast to double.
#
# A factor is rejected rather than coerced: as.numeric() on a factor returns
# its level codes (1, 2, 3, ...), not the values the labels spell, so a count
# stored as a factor would be analyzed as a different count with no warning.
# A character vector is rejected for the same reason: the method needs the
# outcome's numbers, and the conversion is the user's decision. A logical
# outcome is accepted, since TRUE and FALSE coerce to the 1 and 0 a binary
# outcome needs. Any other vector that is not numeric (complex numbers, a
# Date, a difftime) is rejected too, since as.numeric() would drop an
# imaginary part or read a date as a count of days. A matrix with more than
# one column is rejected here, so the message names the outcome instead of
# blaming a row-count mismatch on `X`.
#' @keywords internal
#' @noRd
.as_numeric_vector <- function(x, name) {
  if (is.null(x) || length(x) == 0L)
    stop(sprintf("`%s` must be a non-empty numeric vector (the outcome).",
                 name), call. = FALSE)
  if (is.matrix(x) && ncol(x) != 1L)
    stop(sprintf(paste0("`%s` must be a single outcome: a numeric vector or a ",
                        "one-column matrix, not a matrix with %d columns."),
                 name, ncol(x)), call. = FALSE)
  if (is.matrix(x)) x <- as.vector(x)
  if (is.factor(x))
    stop(sprintf(paste0("`%s` is a factor, and a factor's numeric codes are ",
                        "not its values. Convert it first, for example with ",
                        "as.numeric(as.character(%s)) for a count or ",
                        "continuous outcome, or code a binary outcome as 0/1."),
                 name, name), call. = FALSE)
  if (is.character(x))
    stop(sprintf(paste0("`%s` is a character vector; supply the outcome as ",
                        "numbers, for example with as.numeric(%s)."),
                 name, name), call. = FALSE)
  if (!is.atomic(x) || is.list(x) || !(is.numeric(x) || is.logical(x)))
    stop(sprintf("`%s` must be a numeric vector.", name), call. = FALSE)
  as.numeric(x)
}

# A logical flag must be a single TRUE or FALSE. isTRUE() alone would read
# "yes", NA, 1, or c(TRUE, FALSE) as FALSE without a word, so a user who
# asked for the behavior would silently not get it.
#' @keywords internal
#' @noRd
.validate_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  invisible(TRUE)
}

# match.arg() with an error that names the argument. match.arg()'s own
# message says "'arg' should be one of ...", which does not tell the user
# which argument was wrong. Valid input behaves exactly as under match.arg()
# (the first choice when the argument was left at its default, and unique
# partial matching), because the call is delegated to it; only the error is
# rewritten. `choices` defaults to the calling function's formal default for
# `name`, as in match.arg(), so the list of choices is written once. With
# several.ok = TRUE every entry must identify a choice: match.arg() would
# drop an entry that matches nothing (a misspelling, an ambiguous
# abbreviation, or NA) and return the rest, so such an entry stops here.
#' @keywords internal
#' @noRd
.poemed_match_arg <- function(arg, choices, name = deparse(substitute(arg)),
                              several.ok = FALSE) {
  force(name)
  if (missing(choices)) {
    caller <- sys.function(sys.parent())
    choices <- eval(formals(caller)[[name]], envir = parent.frame())
  }
  refuse <- function() {
    got <- paste(deparse(arg, width.cutoff = 60L, nlines = 1L),
                 collapse = "")
    stop(sprintf("`%s` should be one of %s%s; you supplied %s.", name,
                 paste0("\"", choices, "\"", collapse = ", "),
                 if (several.ok) ", or several of them" else
                   " (a single value)",
                 got), call. = FALSE)
  }
  if (several.ok && is.character(arg) &&
      any(pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE) == 0L))
    refuse()
  tryCatch(match.arg(arg, choices, several.ok = several.ok),
           error = function(e) refuse())
}

# The tuning-parameter grid must be a non-empty vector of strictly positive,
# finite values: each entry is a candidate lambda for the penalized fit, and a
# nonpositive or non-finite lambda is not a valid penalty level.
#' @keywords internal
#' @noRd
.validate_lambda_grid <- function(x, name) {
  if (!is.numeric(x) || length(x) == 0L || any(!is.finite(x)) || any(x <= 0))
    stop(sprintf("`%s` must be a non-empty numeric vector with all elements > 0.",
                 name), call. = FALSE)
  invisible(TRUE)
}

# The tuning grid a worker searches is the `lambda_grid` argument itself; its
# default, seq(0.05, 10, length.out = 100), sits in every signature so it can
# be edited in place. A continuous fit searches the same grid for the full
# model and for the reduced model of the benchmark Wald test, matching the
# reference implementation's default `lamb_grid0 = lamb_grid`. The grid is
# validated by .validate_lambda_grid() above.

# The target error rate for the mediator-screening step (FWER for Bonferroni,
# FDR for BH / BY) must be a single number strictly inside (0, 1): 0 would
# screen out every mediator and 1 would screen out none, and neither endpoint
# is a meaningful target rate.
#' @keywords internal
#' @noRd
.validate_error_level <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1)
    stop("`error_level` must be a single number in (0, 1).", call. = FALSE)
  invisible(TRUE)
}

# Shared structural validation of (X, Y, M, Z): conforming row counts, a
# full-rank unpenalized design [X, Z] (the Wald inference inverts cross
# products of this block, so a collinear column makes the variance estimate
# singular), and no constant mediator columns (a zero-variance mediator cannot
# carry a signal and breaks standardization). Returns the coerced pieces plus
# n, q, and s = ncol(Z). Y-type checks (continuous vs 0/1 vs count) are left to
# the family worker, since only it knows the outcome model.
#' @keywords internal
#' @noRd
.pe_validate_inputs <- function(X, Y, M, Z, drop_constant = FALSE) {
  .validate_flag(drop_constant, "drop_constant")
  # X and M are required. .as_numeric_matrix() passes NULL through (so an
  # absent Z stays absent), so a missing exposure or mediator matrix is caught
  # here, before a row count is read from it.
  if (is.null(X))
    stop("`X` must be supplied: a numeric vector or matrix of exposures.",
         call. = FALSE)
  if (is.null(M))
    stop("`M` must be supplied: a numeric matrix of at least two candidate ",
         "mediators.", call. = FALSE)
  X <- .as_numeric_matrix(X, "X")
  M <- .as_numeric_matrix(M, "M")
  Y <- .as_numeric_vector(Y, "Y")
  n <- length(Y)
  if (ncol(X) == 0L)
    stop("`X` has no columns; supply at least one exposure.", call. = FALSE)
  # Column names, when the user supplied them (a named matrix or a
  # data.frame), travel to the output: the selected set, the mediator table,
  # and the print footer are then read without a lookup against colnames(M),
  # and the per-exposure rows are named after the exposure columns. Captured
  # before any constant column is dropped, so they index the original `M`.
  mediator_names <- .pe_usable_names(M)
  exposure_names <- .pe_usable_names(X)
  # An NA column name is not a usable label (the output falls back to
  # positions), and left on the matrix it breaks summary() of the
  # mediator regression inside the fit, so drop the names from `M`. The
  # names never enter a computed value, so no number changes.
  if (anyNA(colnames(M))) colnames(M) <- NULL

  if (nrow(X) != n)
    stop("`X` must have the same number of rows as the length of `Y`.",
         call. = FALSE)
  if (nrow(M) != n)
    stop("`M` must have the same number of rows as the length of `Y`.",
         call. = FALSE)

  if (!is.null(Z)) {
    Z <- .as_numeric_matrix(Z, "Z")
    if (nrow(Z) != n)
      stop("`Z` must have the same number of rows as the length of `Y`.",
           call. = FALSE)
  }

  # Every value must be finite: a missing, NaN, or infinite entry would
  # otherwise surface as an opaque error deep inside glmnet or the coordinate
  # descent, or silently poison the standardization. The method has no
  # missing-data mechanism; complete cases are the user's decision.
  for (nm in c("X", "Y", "M", "Z")) {
    v <- get(nm)
    if (!is.null(v) && !all(is.finite(v))) {
      bad <- sum(!is.finite(v))
      stop(sprintf(paste0("`%s` contains %d missing, NaN, or infinite value%s. ",
                          "POEMED has no missing-data mechanism; supply complete ",
                          "cases (for example with stats::complete.cases())."),
                   nm, bad, if (bad == 1L) "" else "s"), call. = FALSE)
    }
  }

  # At least two candidate mediators: the penalized selection, the screen,
  # and the multiplicity adjustment are all defined over a set of mediators.
  if (ncol(M) < 2L)
    stop("`M` must have at least two columns (candidate mediators); a single ",
         "mediator calls for an ordinary mediation model.", call. = FALSE)

  # The exposures and confounders must vary: a constant column cannot be
  # standardized and would make the unpenalized design singular after scaling.
  for (nm in c("X", "Z")) {
    v <- get(nm)
    if (!is.null(v)) {
      const <- which(apply(v, 2L, function(col) length(unique(col)) < 2L))
      if (length(const))
        stop(sprintf("`%s` has constant column%s: %s. Remove %s; an intercept is fitted where the model needs one.",
                     nm, if (length(const) == 1L) "" else "s",
                     paste(const, collapse = ", "),
                     if (length(const) == 1L) "it" else "them"), call. = FALSE)
    }
  }

  # The unpenalized design [X, Z] must be full column rank. The mediation
  # inference inverts cross products formed from these columns; a rank
  # deficiency (a duplicated or collinear covariate) would make those
  # inversions fail or return garbage, so we stop early with a clear message.
  W <- if (is.null(Z)) X else cbind(X, Z)
  if (qr(W)$rank < ncol(W))
    stop("The design matrix formed by `X` and `Z` is singular (collinear ",
         "columns). Please remove redundant variables.", call. = FALSE)

  # A mediator with zero variance carries no information and cannot be
  # standardized (its standard deviation is zero). By default this is an error;
  # with drop_constant = TRUE the constant columns are dropped (with a warning)
  # and `kept` records the surviving columns' positions in the original `M`, so
  # the caller can map reported mediator indices back to the user's columns.
  m_var <- apply(M, 2L, stats::var)
  kept <- NULL
  if (any(m_var == 0)) {
    const_cols <- which(m_var == 0)
    if (!isTRUE(drop_constant))
      stop("The mediator matrix `M` contains constant (zero-variance) ",
           "columns: ", paste(const_cols, collapse = ", "),
           ". Remove them, or pass `drop_constant = TRUE` to drop them ",
           "automatically.", call. = FALSE)
    if (length(const_cols) >= ncol(M))
      stop("Every column of the mediator matrix `M` is constant ",
           "(zero variance); there is nothing to test.", call. = FALSE)
    kept <- which(m_var != 0)
    # The two-mediator guard above ran before the drop, so repeat it: one
    # surviving mediator is no longer a high-dimensional mediation problem,
    # and glmnet refuses a one-column matrix. Two columns is only the floor
    # this guard enforces, not a size at which every fit is known to work:
    # with two or three candidate mediators the linear and Poisson fits can
    # still stop inside glmnet or inside the worker, with an error that
    # names no argument.
    if (length(kept) < 2L)
      stop("`M` has only one non-constant column once the constant ",
           "column", if (length(const_cols) > 1L) "s" else "", " (",
           paste(const_cols, collapse = ", "), ") ",
           if (length(const_cols) > 1L) "are" else "is",
           " dropped by `drop_constant = TRUE`; at least two candidate ",
           "mediators are needed, and a single mediator calls for an ",
           "ordinary mediation model.", call. = FALSE)
    warning(length(const_cols), " constant (zero-variance) mediator column",
            if (length(const_cols) > 1L) "s" else "", " dropped: ",
            paste(const_cols, collapse = ", "), ".", call. = FALSE)
    M <- M[, kept, drop = FALSE]
  }

  list(X = X, Y = Y, M = M, Z = Z, n = n, q = ncol(X),
       s = if (is.null(Z)) 0L else ncol(Z), kept = kept,
       mediator_names = mediator_names, exposure_names = exposure_names)
}

# Column names usable as labels: present, non-missing, non-empty, and
# unique. Anything else (no names, an NA name, a blank, a duplicate) falls
# back to column positions. The NA test comes first because nzchar(NA) is
# TRUE, so an NA name would otherwise pass as a label.
#' @keywords internal
#' @noRd
.pe_usable_names <- function(x) {
  nm <- colnames(x)
  if (is.null(nm) || anyNA(nm) || !all(nzchar(nm)) || anyDuplicated(nm))
    return(NULL)
  nm
}

# Validation of a simulation design (simulate_mediation_data()). Each
# argument is checked by name so a bad value stops with a clear message
# instead of surfacing from the covariance construction or rnorm(). Not
# exported.
#' @keywords internal
#' @noRd
.validate_sim_design <- function(n, p, q, d, c1, c2, rho, sigma_y) {
  whole <- function(v, nm, at_least) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < at_least ||
        v != round(v))
      stop(sprintf("`%s` must be a single whole number of at least %d.",
                   nm, at_least), call. = FALSE)
  }
  scalar <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v))
      stop(sprintf("`%s` must be a single finite number.", nm), call. = FALSE)
  }
  whole(n, "n", 2L); whole(p, "p", 2L); whole(q, "q", 1L); whole(d, "d", 0L)
  scalar(c1, "c1"); scalar(c2, "c2"); scalar(rho, "rho"); scalar(sigma_y, "sigma_y")
  if (abs(rho) >= 1)
    stop("`rho` must lie strictly between -1 and 1.", call. = FALSE)
  if (sigma_y <= 0)
    stop("`sigma_y` must be positive.", call. = FALSE)
  invisible(TRUE)
}

# Standardize covariates and mediators (and center the continuous outcome)
# before fitting, following the preprocessing described in Yu and Kelley (in press): the
# penalty treats every mediator coefficient on a common scale, which requires
# the columns to share a scale. X, M, and Z are standardized to mean 0 and
# unit variance. For a continuous (Gaussian) outcome Y is centered so the
# intercept can be dropped from the linear model; for binary and count
# outcomes Y is a 0/1 or count response and is left untouched, since centering
# would destroy its meaning under the link.
#' @keywords internal
#' @noRd
.pe_scale_data <- function(X, Y, M, Z, center_y) {
  X <- base::scale(X)
  M <- base::scale(M)
  if (!is.null(Z)) Z <- base::scale(Z)
  if (center_y) Y <- Y - mean(Y)
  list(X = X, Y = Y, M = M, Z = Z)
}

# The global test rejects when the test statistic exceeds the upper-alpha
# quantile of a chi-square with q degrees of freedom (q = number of
# exposures), so the p-value is the chi-square upper-tail probability. Both the
# benchmark Wald statistic S_n and the power-enhanced statistic M_PE use this
# same reference distribution under the global null (Theorems on the limiting
# null distribution). lower.tail = FALSE is used rather than 1 - pchisq() so
# the far upper tail (p-values near machine zero, which the PE statistic
# routinely produces) keeps its precision instead of canceling to 0.
#' @keywords internal
#' @noRd
.pe_chisq_pvalue <- function(stat, df) {
  stats::pchisq(stat, df = df, lower.tail = FALSE)
}

# Optionally-parallel lapply used by the Monte Carlo and multi-group routines.
# On Unix with cores > 1 it forks via parallel::mclapply (parallel ships with
# R, so it adds no installation, and it is declared in Imports because it is
# called with ::); otherwise it is an ordinary lapply. Callers
# that need reproducibility under forking set RNGkind("L'Ecuyer-CMRG") and a
# seed first, which gives parallel-safe streams. Not exported.
#' @keywords internal
#' @noRd
.pe_lapply <- function(X, FUN, cores = 1L) {
  if (cores > 1L && .Platform$OS.type == "unix")
    parallel::mclapply(X, FUN, mc.cores = cores)
  else
    lapply(X, FUN)
}
