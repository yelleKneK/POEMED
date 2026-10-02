#' Tidy Printing for POEMED Result Tables
#'
#' Every estimation and testing function in \pkg{POEMED} returns a tidy
#' \code{data.frame} with a \code{term} column and one or more numeric
#' columns. A single numeric column routinely holds quantities on very
#' different scales: a whole-number count of mediators next to a
#' chi-square statistic in the hundreds next to a p-value of
#' \eqn{10^{-12}}{10^-12}. The base \code{\link[base]{print.data.frame}} method
#' formats a whole column with one common format, which forces either a
#' wall of trailing zeros or a slide into scientific notation. The
#' \code{poemed_tbl} class supplies \code{print} and \code{format} methods
#' that format each value on its own terms: whole numbers (counts,
#' sample sizes, mediator indices) print without a decimal part, other
#' values print to a few significant figures, and p-values print to a
#' fixed number of decimal places with a \dQuote{< 0.0001} floor.
#'
#' The stored numeric values are never rounded. Only the display
#' changes, so any arithmetic you do on the returned object (a
#' difference of two statistics, a further calculation) uses full
#' precision. To see more digits, raise \code{digits} or read the column
#' directly with \code{x$value}.
#'
#' This mirrors the display convention of the \pkg{DMAR} package, whose
#' house style \pkg{POEMED} follows.
#'
#' @section Subsetting and Combining:
#' A \code{poemed_tbl} carries its display information as attributes: which
#' rows of the \code{value} column hold p-values, the outcome model, the
#' selected mediators, and the tuning record that the print footer
#' reports. Subsetting with \code{[} or \code{\link[base]{subset}()}
#' keeps all of them, whether you select rows, columns, or both, so a
#' subset still prints its p-values with the floor and names the fit it
#' came from. (Base R would drop them whenever columns are selected.)
#'
#' Combining tables with \code{rbind()} keeps these attributes only when
#' every argument is a \code{poemed_tbl} and all of them carry the same
#' ones, as when you split one table and put it back together. Tables from
#' different fits describe different models, so their combination keeps
#' only the list of p-value rows and prints no single-fit footer. A row
#' whose \code{term} begins with \code{pval_} prints as a p-value even
#' when that list is missing.
#'
#' @param x A \code{poemed_tbl} object (a tidy \code{data.frame} returned
#'   by a POEMED function).
#' @param digits Number of significant figures for non-integer values, a
#'   single whole number from 1 to 22 (the range base R's
#'   \code{\link[base]{format}()} accepts). Defaults to
#'   \code{getOption("poemed.digits", 4L)}, so
#'   \code{options(poemed.digits = 6)} changes it for the session.
#' @param digits_p Number of decimal places for p-values, a single whole
#'   number from 1 to 15 (about the number of decimal digits a double holds
#'   reliably). Defaults to 4. A p-value below \code{10^(-digits_p)} prints as
#'   \dQuote{< } followed by that bound, so \dQuote{< 0.0001} at the
#'   default and \dQuote{< 0.01} with \code{digits_p = 2}.
#' @param ... For \code{print}, further arguments passed to
#'   \code{\link[base]{print.data.frame}}; \code{row.names} and
#'   \code{right} default to \code{FALSE} for the tidy look. For
#'   \code{format}, ignored. For \code{[}, the row and column indices and
#'   \code{drop}, as for a data frame. For \code{rbind}, the tables to
#'   combine and any options of \code{\link[base]{rbind.data.frame}}, such
#'   as \code{make.row.names}.
#' @param deparse.level Passed to \code{\link[base]{rbind}()}.
#'
#' @return \code{print.poemed_tbl} returns \code{x} invisibly.
#'   \code{format.poemed_tbl} returns a \code{data.frame} whose numeric
#'   columns have been formatted to character for display. The \code{[}
#'   method returns a \code{poemed_tbl} with the attributes of \code{x}
#'   (or a vector, when a single column is extracted), and
#'   \code{rbind.poemed_tbl} returns a \code{poemed_tbl} whose attributes
#'   follow the rule under Subsetting and Combining.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @examples
#' set.seed(113)
#' d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
#'                              outcome = "continuous")
#' fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
#' # HBIC chose the smallest value of the default grid here, so the fit
#' # warned; see the lambda_grid argument of ?pe_mediation.
#' fit                       # rounded for reading, with the p-value floor
#' fit$value[fit$term == "pval_pe"]   # full precision underneath
#'
#' # A subset keeps the p-value format and the footer of its fit.
#' subset(fit, term %in% c("stat_pe", "pval_pe"))
#'
#' @name poemed_tbl
NULL

# Element-wise display formatter shared by format.poemed_tbl / print.poemed_tbl.
# A finite value that equals its own rounding (a whole number, such as a
# count of mediators or a sample size) prints with no decimal part; any
# other finite value prints to `digits` significant figures, letting base R
# choose fixed versus scientific notation per value (so a chi-square of 412.7
# reads plainly while a tiny variance reads as 1.2e-05). The whole-number
# branch is capped at 1e15 because beyond that doubles cannot represent
# integers exactly and "x == round(x)" stops being meaningful. Not exported.
.format_poemed_value <- function(v, digits = 4L) {
  vapply(v, function(x) {
    if (is.na(x)) return("NA")
    if (is.finite(x) && x == round(x) && abs(x) < 1e15)
      return(format(x, scientific = FALSE, trim = TRUE))
    format(signif(x, digits), digits = digits, trim = TRUE)
  }, character(1L))
}

# Fixed-decimal p-value formatter. p-values read best at a fixed number of
# decimal places (the convention is four), never in scientific notation. A
# value below the smallest representable magnitude (10^(-digits_p)) prints as
# "< 0.0001" rather than rounding to "0.0000", which would misread as an exact
# zero. The PE tests can return p-values near machine zero (the chi-square
# upper tail at a statistic in the hundreds), so the floor matters here. Not
# exported.
.format_poemed_pvalue <- function(p, digits_p = 4L) {
  thresh <- 10^(-digits_p)
  floor_label <- paste0("< ", formatC(thresh, format = "f", digits = digits_p))
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < thresh) return(floor_label)
    formatC(x, format = "f", digits = digits_p)
  }, character(1L))
}

# A digits argument must be a single whole number from 1 to `upper`; anything
# else would make format() or formatC() error obscurely or print nothing.
# The upper bounds are 22 for `digits`, the largest value base R's format()
# accepts (23 fails inside prettyNum() with a message that names no
# argument), and 15 for `digits_p`, about the number of decimal digits a
# double holds reliably (further decimals print binary-expansion noise).
# `option` names the global option the value came from when the caller left
# the argument at its default, so the message can point the user to the
# setting that needs changing. Not exported.
.validate_digits <- function(x, name, upper = 22L, option = NULL) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 ||
      x > upper || x != round(x)) {
    fix <- if (is.null(option)) sprintf("supply, for example, `%s = 4`", name)
           else sprintf(paste0("it defaults to getOption(\"%s\"), so reset ",
                               "that option with, for example, ",
                               "options(%s = 4)"), option, option)
    stop(sprintf("`%s` must be a single whole number from 1 to %d; %s.",
                 name, upper, fix), call. = FALSE)
  }
  invisible(TRUE)
}

# The attributes of a poemed_tbl that describe its content (everything but
# the data frame's own names, row names, and class). Not exported.
.poemed_tbl_meta <- function(x) {
  a <- attributes(x)
  a[setdiff(names(a), c("names", "row.names", "class"))]
}

# Tag a tidy result table as a poemed_tbl so it prints via print.poemed_tbl. Every
# tidy-returning POEMED function routes its result through this helper just
# before returning it. `p_terms`, when supplied, names the rows of the long
# `value` column that hold p-values (so the format method gives them fixed
# decimals); a wide table with a literal `p_value` column needs no p_terms, as
# the format method finds that column by name. The poemed_tbl class is inserted
# just before `data.frame`, after any leading tidy/glance subclass, so the
# subclass keeps dispatching its own methods while print falls through here.
# Idempotent. Not exported.
.as_poemed_tbl <- function(out, p_terms = NULL) {
  if (!is.null(p_terms)) attr(out, "p_terms") <- p_terms
  cls <- class(out)
  if (!("poemed_tbl" %in% cls))
    class(out) <- c(setdiff(cls, "data.frame"), "poemed_tbl", "data.frame")
  out
}

#' @rdname poemed_tbl
#' @export
format.poemed_tbl <- function(x, digits = getOption("poemed.digits", 4L),
                            digits_p = 4L, ...) {
  .validate_digits(digits, "digits", 22L,
                   option = if (missing(digits)) "poemed.digits")
  .validate_digits(digits_p, "digits_p", 15L)
  raw <- x
  class(raw) <- "data.frame"
  out <- raw
  num <- vapply(out, is.numeric, logical(1L))
  out[num] <- lapply(out[num], .format_poemed_value, digits = digits)

  # The "primary" numeric column carries the per-term quantities in a long
  # table: `value`. The p_terms attribute names rows of this column whose
  # quantity is a p-value and so should print to fixed decimals.
  primary <- intersect("value", names(out)[num])
  primary <- if (length(primary)) primary[1L] else NA_character_

  # p-value columns in a wide table are detected by name: p_value, p.value,
  # pval_* (the fit and WHO tables), and screen_p (the mediator table).
  pcols <- grep("^(p_value|p\\.value|pval(_.*)?|screen_p)$", names(out)[num],
                value = TRUE)
  for (nm in pcols)
    out[[nm]] <- .format_poemed_pvalue(raw[[nm]], digits_p = digits_p)

  # Named rows of the long-format `value` column are detected via p_terms.
  # A long table that has lost that attribute (say, after a merge) still has
  # its pval_* rows recognized by name, since in POEMED a term of that name
  # always holds a p-value.
  p_terms <- attr(x, "p_terms")
  if (is.null(p_terms) && "term" %in% names(raw))
    p_terms <- grep("^pval(_|$)", as.character(raw$term), value = TRUE)
  if (length(p_terms) && "term" %in% names(out) && !is.na(primary)) {
    idx <- raw$term %in% p_terms
    out[[primary]][idx] <- .format_poemed_pvalue(raw[[primary]][idx],
                                               digits_p = digits_p)
  }
  out
}

#' @rdname poemed_tbl
#' @export
print.poemed_tbl <- function(x, digits = getOption("poemed.digits", 4L),
                           digits_p = 4L, ...) {
  # Validate here as well as in format() so that a bad value set through
  # options(poemed.digits) is reported as coming from that option.
  .validate_digits(digits, "digits", 22L,
                   option = if (missing(digits)) "poemed.digits")
  .validate_digits(digits_p, "digits_p", 15L)
  disp <- format(x, digits = digits, digits_p = digits_p)
  # row.names / right travel through `...` so a caller can override them, but
  # default to the tidy look (no row numbers, left-aligned).
  dots <- list(...)
  if (is.null(dots[["row.names"]])) dots[["row.names"]] <- FALSE
  if (is.null(dots[["right"]]))     dots[["right"]]     <- FALSE
  do.call(print.data.frame, c(list(disp), dots))

  # A short, human-readable footer naming the outcome model and the
  # mediators the screen selected. The numbers are already in the table; the
  # footer just makes the headline reading effortless.
  # A fit names its outcome model; a study or WHO table names its outcome.
  outcome <- attr(x, "outcome")
  if (!is.null(outcome))
    cat(sprintf("\n%s: %s\n", if ("term" %in% names(x)) "Outcome model" else "Outcome",
                outcome))
  sel <- attr(x, "selected_mediators")
  if (!is.null(sel)) {
    if (length(sel) == 0L) {
      cat("Selected mediators: none\n")
    } else {
      # By column name when `M` had names, by column position otherwise.
      nm <- attr(x, "mediator_names")
      lab <- if (is.null(nm)) sel else nm[sel]
      cat(sprintf("Selected mediators (%d): %s\n",
                  length(sel), paste(lab, collapse = ", ")))
    }
  }
  # The tuning record: which lambda HBIC chose from which grid, and whether it
  # sat at either end of the grid (HBIC is minimized over the grid alone, so
  # its minimum may lie beyond that end; the fit warned when that happened).
  # The footer reports the full model's choice; for a continuous outcome the
  # reduced model's choice on the same grid, `lambda_selected_reduced`, is in
  # the attribute.
  tun <- attr(x, "tuning")
  if (!is.null(tun)) {
    g <- tun$lambda_grid
    cat(sprintf("Tuning parameter (HBIC): lambda = %s from %d values in [%s, %s]%s\n",
                format(signif(tun$lambda_selected, 3)), length(g),
                format(signif(min(g), 3)), format(signif(max(g), 3)),
                if (isTRUE(tun$at_lower_end)) " (the grid's lower end)"
                else if (isTRUE(tun$at_upper_end)) " (the grid's upper end)"
                else ""))
  }
  # An empty fit means the HBIC-selected penalized fit kept no mediator;
  # smaller grid values may have kept some and lost on HBIC, so the footer
  # names the selected fit rather than claiming every lambda was empty.
  if (isTRUE(attr(x, "empty_fit")))
    cat(paste0("Penalized selection: empty at the HBIC-selected lambda; ",
               "no test was computed.\n"))
  invisible(x)
}

#' @rdname poemed_tbl
#' @export
`[.poemed_tbl` <- function(x, ...) {
  out <- NextMethod()
  # A single extracted column is a plain vector; return it as base R does.
  if (!is.data.frame(out)) return(out)
  # Base R keeps the attributes on a row subset but drops them when columns
  # are selected; restore any that were dropped so the subset keeps its
  # p-value rows, footer, and outcome label.
  meta <- .poemed_tbl_meta(x)
  for (nm in names(meta))
    if (is.null(attr(out, nm))) attr(out, nm) <- meta[[nm]]
  out
}

#' @rdname poemed_tbl
#' @method rbind poemed_tbl
#' @export
rbind.poemed_tbl <- function(..., deparse.level = 1) {
  args <- list(...)
  nms <- names(args)
  if (is.null(nms)) nms <- rep("", length(args))
  # Named options of rbind.data.frame() travel through `...` with the tables.
  is_opt <- nms %in% c("make.row.names", "stringsAsFactors", "factor.exclude")
  parts <- args[!is_opt]
  tbls <- Filter(function(a) inherits(a, "poemed_tbl"), parts)
  meta <- lapply(tbls, function(a) {
    m <- .poemed_tbl_meta(a)
    m[order(names(m))]
  })
  # Strip every table to a plain data frame so rbind.data.frame() does not
  # copy the first argument's attributes onto the whole result.
  plain <- lapply(parts, function(a) {
    if (!inherits(a, "poemed_tbl")) return(a)
    attributes(a) <- attributes(a)[c("names", "row.names")]
    class(a) <- "data.frame"
    a
  })
  out <- do.call(rbind, c(plain, args[is_opt],
                          list(deparse.level = deparse.level)))
  # Keep the attributes only when every piece is a poemed_tbl carrying the
  # same ones (a table split and put back together). Otherwise the pieces
  # describe different fits, and no single footer is true of all the rows;
  # keep just the union of the p-value rows so they still print with the
  # floor.
  same <- length(tbls) > 0L && length(tbls) == length(parts) &&
    all(vapply(meta, identical, logical(1L), meta[[1L]]))
  if (same) {
    for (nm in names(meta[[1L]])) attr(out, nm) <- meta[[1L]][[nm]]
    p_terms <- attr(out, "p_terms")
  } else {
    p_terms <- unique(unlist(lapply(tbls, attr, "p_terms")))
  }
  .as_poemed_tbl(out, p_terms = p_terms)
}
