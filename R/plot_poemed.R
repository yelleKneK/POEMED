#' Plot a POEMED Result
#'
#' A base-graphics plot for the two kinds of table POEMED returns. For a
#' power curve from [pe_power_curve()] it draws the empirical rejection
#' rate of the benchmark and power-enhanced tests against the
#' signal-strength grid, with the nominal level marked. For a single
#' [pe_mediation()] fit it draws the per-mediator screening evidence
#' (\eqn{-\log_{10}}{-log10} of the screening p-value) for each candidate
#' mediator the penalized fit kept, with the selected ones highlighted, so it is clear
#' which mediators drove the result.
#'
#' A screening p-value is stored as exactly 0 only when it underflows the
#' smallest double, which takes both path statistics at about 38 or more
#' in absolute value (the table prints it as `< 0.0001`, like any p-value
#' below that floor). Its
#' \eqn{-\log_{10}}{-log10} is infinite, so the per-mediator view draws such a
#' point at a cap instead. The cap is the larger of \eqn{-\log_{10}}{-log10} of
#' the machine epsilon (about 15.65) and the tallest finite point. A capped
#' point is drawn as a triangle, a dotted line marks the cap, and the
#' legend says so. The true height of a capped point is at least the cap.
#' Only the drawing changes; the stored `screen_p` values (see
#' [pe_mediators()]) are untouched.
#'
#' @param x A `poemed_tbl` from [pe_power_curve()] or [pe_mediation()].
#' @param ... Further graphical parameters passed to the underlying
#'   [plot.default()] call. They override the defaults, which are `xlab`,
#'   `ylab`, `pch`, `col`, and `ylim` (from 0 to 1.3 times the tallest
#'   point, leaving room for the legend) in the per-mediator view and
#'   `xlab`, `ylab`, `pch`, `type`, and `ylim` in the power-curve view, so
#'   `plot(fit, xlab = "mediator", col = "red")` works. A supplied `pch`,
#'   `col`, or `lty` is carried into the legend. In the per-mediator view
#'   the x axis is labeled with the mediators' names or positions; passing
#'   `xaxt` replaces that axis with the one `xaxt` asks for.
#'
#' @return `x`, invisibly. Called for the plot it draws. A table that is
#'   neither a fit with a per-mediator table nor a power curve draws
#'   nothing and gives a message.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_power_curve()], [pe_mediation()], [pe_mediators()].
#'
#' @family mediation simulation
#'
#' @examples
#' set.seed(113)
#' d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
#'                              pattern = "contrasting", c1 = 1)
#' fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
#' # HBIC chose the smallest value of the default grid here, so the fit
#' # warned; see the lambda_grid argument of ?pe_mediation.
#' plot(fit)
#' plot(fit, main = "Screening evidence", xlab = "candidate mediator", pch = 15)
#'
#' @method plot poemed_tbl
#' @export
plot.poemed_tbl <- function(x, ...) {
  dots <- list(...)
  if (all(c("c1", "rejection_hdmm", "rejection_pe") %in% names(x))) {
    # Power-curve view: two rejection-rate curves against the signal grid.
    # The defaults are merged with the user's graphical parameters (theirs
    # win) and passed in one call, so setting xlab or ylim does not collide
    # with a default of the same name.
    alpha_level <- attr(x, "alpha_level")
    args <- utils::modifyList(
      list(x = x$c1, y = x$rejection_pe, type = "b", pch = 19,
           ylim = c(0, 1), xlab = "signal strength (c1)",
           ylab = "empirical rejection rate"),
      dots)
    do.call(graphics::plot, args)
    graphics::lines(x$c1, x$rejection_hdmm, type = "b", pch = 1, lty = 2)
    if (!is.null(alpha_level)) graphics::abline(h = alpha_level, col = "grey60", lty = 3)
    fg <- graphics::par("col")
    graphics::legend("right", c("power-enhanced", "benchmark Wald"),
                     pch = c(if (is.null(args$pch)) 19 else args$pch[1L], 1),
                     lty = c(if (is.null(args$lty)) 1 else args$lty[1L], 2),
                     col = c(if (is.null(args$col)) fg else args$col[1L], fg),
                     bty = "n")
    return(invisible(x))
  }

  tab <- attr(x, "mediator_table")
  if (!is.null(tab) && nrow(tab) > 0L) {
    # Per-mediator view: screening evidence, selected mediators highlighted.
    # A p-value stored as exactly 0 is drawn at a marked cap rather than at
    # an infinite (or, clamped to the smallest double, a 307-decade) height.
    ev <- .screen_evidence(tab$screen_p)
    k <- nrow(tab)
    # The y axis starts at 0 (a p-value of 1) and leaves headroom above the
    # tallest point, where the legend sits, so the legend does not cover the
    # strongest mediators (capped ones all share the top height).
    top <- suppressWarnings(max(ev$height, na.rm = TRUE))
    ylim <- c(0, if (is.finite(top) && top > 0) 1.3 * top else 1)
    args <- utils::modifyList(
      list(x = seq_len(k), y = ev$height,
           pch = ifelse(ev$capped, 17, 19),
           col = ifelse(tab$selected, "black", "grey70"), xaxt = "n",
           ylim = ylim, xlab = "selected candidate mediator",
           ylab = expression(-log[10](screening~p))),
      dots)
    do.call(graphics::plot, args)
    if (is.null(dots$xaxt))
      graphics::axis(1, at = seq_len(k),
                     labels = if ("name" %in% names(tab)) tab$name else tab$mediator)
    if (any(ev$capped)) graphics::abline(h = ev$cap, col = "grey60", lty = 3)
    # The legend takes each entry's symbol and color from a point of its
    # kind (falling back to the default when there is none), so it stays
    # true when the user supplies pch or col. A capped point has a screening
    # p-value of 0, which clears any threshold, so it is always selected.
    pch_all <- rep_len(if (is.null(args$pch)) graphics::par("pch") else args$pch, k)
    col_all <- rep_len(if (is.null(args$col)) graphics::par("col") else args$col, k)
    from <- function(v, which, default) if (any(which)) v[which][1L] else default
    def_pch <- if (is.null(dots$pch)) 19 else pch_all[1L]
    def_col <- function(d) if (is.null(dots$col)) d else col_all[1L]
    act <- tab$selected
    cap <- ev$capped
    leg     <- c("selected", "screened out")
    leg_pch <- c(from(pch_all, act & !cap, def_pch),
                 from(pch_all, !act & !cap, def_pch))
    leg_col <- c(from(col_all, act, def_col("black")),
                 from(col_all, !act, def_col("grey70")))
    if (any(cap)) {
      leg     <- c(leg, "p stored as 0 (drawn at cap)")
      leg_pch <- c(leg_pch, from(pch_all, cap, 17))
      leg_col <- c(leg_col, from(col_all, cap, def_col("black")))
    }
    graphics::legend("topright", leg, pch = leg_pch, col = leg_col, bty = "n")
    return(invisible(x))
  }

  message("No POEMED plot is defined for this table.")
  invisible(x)
}

# Heights for the per-mediator view: -log10 of each screening p-value, with a
# p-value stored as exactly 0 drawn at a cap instead of at infinity. A stored
# 0 is a p-value that underflowed the smallest double: the screen computes
# 2 * pnorm(-|t|) (see pe_component.R), which is nonzero until both path
# statistics pass about 38 in absolute value, so a stored 0 has a true
# height above 323 decades. The cap is -log10(machine epsilon), about 15.65,
# or the tallest finite height if that is larger, so a capped point never
# sits below a point whose p-value was stored as nonzero. Returns the
# heights, which points are capped, and the cap. The stored p-values are not
# modified. Not exported.
.screen_evidence <- function(screen_p) {
  capped <- !is.na(screen_p) & screen_p <= 0
  height <- rep(NA_real_, length(screen_p))
  height[!capped] <- -log10(screen_p[!capped])
  finite <- height[is.finite(height)]
  cap <- max(-log10(.Machine$double.eps), finite)
  height[capped] <- cap
  list(height = height, capped = capped, cap = cap)
}

#' Broom Verbs for POEMED Results
#'
#' `tidy()` and `glance()` methods so a POEMED result composes with the
#' broom ecosystem; the generics are re-exported, so `tidy(fit)` and
#' `glance(fit)` work after `library(POEMED)` alone. For a [pe_mediation()]
#' fit, `tidy()` returns the per-mediator table (see [pe_mediators()]) and
#' `glance()` returns a one-row summary of the two tests and the total
#' indirect effect. For a [pe_power_curve()] table both return the table
#' itself (it is already one row per grid point).
#'
#' @param x A `poemed_tbl` result.
#' @param ... Unused, for generic compatibility.
#'
#' @return A `data.frame`. For a fit, `tidy()` has the columns of
#'   [pe_mediators()] and `glance()` the columns `stat_hdmm`, `pval_hdmm`,
#'   `stat_pe`, `pval_pe`, `total_indirect_effect` (one column per
#'   exposure when there are several, suffixed with the exposure's column
#'   name or, for an unnamed `X`, `_1`, `_2`, ...), `n_selected_mediators`,
#'   and `n_observations`.
#'
#' @examples
#' set.seed(113)
#' d <- simulate_mediation_data(n = 120, p = 40, pattern = "contrasting",
#'                              outcome = "continuous", c1 = 1)
#' fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
#' # HBIC chose the smallest value of the default grid here, so the fit
#' # warned; see the lambda_grid argument of ?pe_mediation.
#' tidy(fit)
#' glance(fit)
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()], [pe_mediators()].
#'
#' @name poemed_broom
#' @importFrom generics tidy
#' @export
tidy.poemed_tbl <- function(x, ...) {
  tab <- attr(x, "mediator_table")
  if (!is.null(tab)) {
    # Plain row numbers: the internal row names (M.1, M.2, ...) come from the
    # screening regression and need not match the `mediator` column.
    tab <- as.data.frame(tab)
    rownames(tab) <- NULL
    return(tab)
  }
  as.data.frame(unclass(x))
}

#' @rdname poemed_broom
#' @importFrom generics glance
#' @export
glance.poemed_tbl <- function(x, ...) {
  if (!"term" %in% names(x)) return(as.data.frame(unclass(x)))
  v <- function(t) {
    hit <- x$value[x$term == t]
    if (length(hit)) hit[1] else NA_real_
  }
  # One total-indirect-effect column per exposure (the rows are suffixed
  # _1, _2, ... when there are several exposures).
  beta_terms <- grep("^total_indirect_effect(_.+)?$", x$term, value = TRUE)
  beta <- stats::setNames(as.list(x$value[match(beta_terms, x$term)]), beta_terms)
  cbind(data.frame(stat_hdmm = v("stat_hdmm"), pval_hdmm = v("pval_hdmm"),
                   stat_pe = v("stat_pe"), pval_pe = v("pval_pe"),
                   stringsAsFactors = FALSE),
        as.data.frame(beta, stringsAsFactors = FALSE),
        data.frame(n_selected_mediators = v("n_selected_mediators"),
                   n_observations = v("n_observations"),
                   stringsAsFactors = FALSE))
}

#' @rdname poemed_broom
#' @export
generics::tidy

#' @rdname poemed_broom
#' @export
generics::glance
