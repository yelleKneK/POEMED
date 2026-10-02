#' Summarize a POEMED Table
#'
#' For a grouped table that compares the benchmark and power-enhanced tests
#' (any `poemed_tbl` with a `group` column, a `pval_hdmm` column, and a
#' power-enhanced p-value column (`pval_pe` or `pval_pe_bonferroni`), such as
#' the output of [WHO_mediation_analysis()]), `summary()` digests it across
#' groups: how many groups each test detects mediation in, and in particular
#' the groups where the power-enhanced test detects mediation that the
#' benchmark misses, together with the most frequently flagged mediators. For
#' any other `poemed_tbl` it falls back to the ordinary data-frame summary.
#'
#' The printed digest follows the display convention of [print.poemed_tbl()]:
#' p-values appear to four decimal places, with a `< 0.0001` floor and never
#' in scientific notation, and the significance level appears as supplied.
#' The returned object keeps the p-values at full precision.
#'
#' @param object A `poemed_tbl`.
#' @param alpha_level Significance level for counting a group as a detection.
#'   Default 0.05.
#' @param ... For a grouped comparison table, ignored. For any other
#'   `poemed_tbl`, passed to [summary.data.frame()] (for example `digits` or
#'   `quantile.type`).
#' @param x A `summary.poemed_tbl` object, for the print method.
#'
#' @return For a grouped comparison table, an object of class
#'   `summary.poemed_tbl` (a list, printed by its own method) with elements
#'   `outcome` (the outcome label, or `NULL` when the table has none),
#'   `full_table` (logical; `TRUE` for the extended layout of
#'   [WHO_mediation_analysis()] with `full_table = TRUE`), `alpha_level`,
#'   `n_groups` (rows in the table), `n_with_data` (groups with a benchmark
#'   p-value), `n_hdmm` and `n_pe` (groups detected by each test),
#'   `pe_only` (a `data.frame` of the groups the PE test detects but the
#'   benchmark does not, with their `n_countries`, both p-values, and the
#'   `selected_mediators` codes), and `top_mediators` (a frequency table of
#'   the mediators selected in the groups the PE test detects). For any other
#'   `poemed_tbl`, the value of [summary.data.frame()]. The print method
#'   returns its argument invisibly.
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [WHO_mediation_analysis()], [pe_mediation()].
#'
#' @family mediation tests
#'
#' @examples
#' # The global and income-group infant-mortality models, on a 25-point
#' # version of the default 100-point grid to keep the example fast (the
#' # global model's grid defaults to the same grid). The coarser grid
#' # selects the same mediators as the default, those of the article's
#' # Table 1, and moves some benchmark p-values slightly.
#' res <- WHO_mediation_analysis("imr", groupings = c("global", "income"),
#'                               lambda_grid = seq(0.1, 2, length.out = 25))
#' summary(res)
#'
#' @exportS3Method summary poemed_tbl
summary.poemed_tbl <- function(object, alpha_level = 0.05, ...) {
  .validate_alpha_level(alpha_level, "alpha_level")
  df <- as.data.frame(object)
  pe_col <- if ("pval_pe" %in% names(df)) "pval_pe"
            else if ("pval_pe_bonferroni" %in% names(df)) "pval_pe_bonferroni"
            else NA_character_
  # Only a grouped benchmark-vs-PE table gets the comparison digest; anything
  # else (a single fit, a power curve) falls back to the standard summary.
  if (!(all(c("group", "pval_hdmm") %in% names(df)) && !is.na(pe_col)))
    return(summary(df, ...))

  pe_p <- df[[pe_col]]
  pe_a <- if ("selected_mediators" %in% names(df)) df$selected_mediators
          else if ("selected_bonferroni" %in% names(df)) df$selected_bonferroni
          else rep(NA_character_, nrow(df))
  n_obs <- if ("n_countries" %in% names(df)) df$n_countries
           else rep(NA_integer_, nrow(df))
  has <- !is.na(df$pval_hdmm)                       # groups with data

  hdmm_sig <- has & df$pval_hdmm <= alpha_level
  pe_sig   <- has & pe_p <= alpha_level
  only     <- pe_sig & !hdmm_sig                    # PE detects, benchmark does not

  pe_only <- data.frame(
    group = df$group[only], n_countries = n_obs[only],
    pval_hdmm = df$pval_hdmm[only], pval_pe = pe_p[only],
    selected_mediators = pe_a[only],
    row.names = NULL, stringsAsFactors = FALSE)

  codes <- unlist(strsplit(pe_a[pe_sig & pe_a != "none" & !is.na(pe_a)], ", "))
  top <- if (length(codes)) sort(table(codes), decreasing = TRUE) else integer(0)

  structure(list(
    outcome = attr(object, "outcome"),
    full_table = isTRUE(attr(object, "full_table")),
    alpha_level = alpha_level,
    n_groups = nrow(df), n_with_data = sum(has),
    n_hdmm = sum(hdmm_sig), n_pe = sum(pe_sig),
    pe_only = pe_only, top_mediators = top),
    class = "summary.poemed_tbl")
}

#' @rdname summary.poemed_tbl
#' @exportS3Method print summary.poemed_tbl
print.summary.poemed_tbl <- function(x, ...) {
  # The header names the outcome when the table carries one.
  lab <- if (is.null(x$outcome)) "" else paste0(": ", toupper(x$outcome))
  cat(sprintf("POEMED comparison%s%s\n", lab,
              if (isTRUE(x$full_table)) " (extended table)" else ""))
  # The level prints as supplied (up to 15 significant digits), so a level
  # such as 0.9999999 is not rounded to a value the validator would reject.
  cat(sprintf("  groups: %d (%d with data), alpha_level = %s\n",
              x$n_groups, x$n_with_data, format(x$alpha_level, digits = 15)))
  cat(sprintf("  detected by benchmark (HDMM): %d   by power-enhanced (PE): %d\n",
              x$n_hdmm, x$n_pe))
  n_only <- nrow(x$pe_only)
  if (n_only == 0L) {
    cat("  PE detects no group the benchmark misses.\n")
  } else {
    cat(sprintf("  PE detects mediation in %d group%s the benchmark misses:\n",
                n_only, if (n_only > 1L) "s" else ""))
    # The p-values print with the same fixed-decimal formatter (and floor)
    # as print.poemed_tbl, so a p-value reads the same in both views; the
    # stored values in x$pe_only keep full precision.
    show <- x$pe_only
    show$pval_hdmm <- .format_poemed_pvalue(show$pval_hdmm)
    show$pval_pe <- .format_poemed_pvalue(show$pval_pe)
    print(show, row.names = FALSE)
  }
  if (length(x$top_mediators)) {
    cat("  most-flagged mediators:\n")
    tm <- x$top_mediators
    for (i in seq_along(tm))
      cat(sprintf("    %-18s %d\n", names(tm)[i], tm[i]))
  }
  invisible(x)
}
