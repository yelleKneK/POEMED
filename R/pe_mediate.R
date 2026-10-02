#' Power-Enhanced Mediation Test From a Data Frame (Formula Interface)
#'
#' A convenience front end to [pe_mediation()] for users who have a data
#' frame rather than ready-made matrices. Give it the outcome and exposure
#' as a formula, name the mediator columns, and (optionally) name
#' confounder columns; it assembles the exposure, outcome, mediator, and
#' confounder matrices (turning factor confounders into indicator
#' variables) and runs the test. This is the gentler entry point; the
#' matrix interface [pe_mediation()] is the workhorse.
#'
#' @param formula A two-sided formula `outcome ~ exposure` naming the
#'   outcome and one or more exposure columns in `data`. A `.` stands for
#'   every other column of `data`, so either name the exposures explicitly
#'   or subtract the mediator and confounder columns (`y ~ . - m1 - m2`); a
#'   `.` that brings a mediator or confounder in as an exposure stops.
#' @param data A `data.frame` containing the outcome, exposure, mediator,
#'   and any confounder columns.
#' @param mediators Character vector of mediator column names in `data`
#'   (the candidate mediators `M`). At least two columns, each listed once,
#'   and none of them the outcome, an exposure, or a confounder.
#' @param confounders Optional character vector of confounder column names.
#'   Factor and character columns are expanded into indicator variables;
#'   numeric columns enter as is. Default `NULL`; an empty vector,
#'   `character(0)`, also means no confounders.
#' @param outcome Outcome type: `"continuous"`, `"binary"`, or `"count"`.
#' @param ... Further arguments passed to [pe_mediation()] (for example
#'   `method`, `error_level`, `report_all_methods`). Among them is
#'   `lambda_grid`, the grid of candidate tuning parameters that HBIC
#'   searches, `seq(0.05, 10, length.out = 100)` by default for every
#'   outcome family: a shorter or narrower grid makes the fit faster, and a
#'   wider or finer one is worth trying when the value chosen is the
#'   smallest or largest in the grid (a non-empty fit warns, and the print
#'   footer, which reports the full model's choice, says so). See the
#'   argument on [pe_mediation()].
#'
#' @details
#' Each column plays one role. The outcome and the exposures come from
#' `formula`, the mediators from `mediators`, and the confounders from
#' `confounders`. A column listed in two roles, or a mediator listed twice,
#' stops with a message naming it: the outcome listed as a mediator, for
#' example, would be "identified" as its own mediator.
#'
#' Missing values follow one complete-case rule across every column the
#' model uses (outcome, exposures, mediators, and confounders). A row with a
#' missing value in any of them is dropped, and a message reports how many
#' rows were dropped and which columns held the missing values. The method
#' has no missing-data mechanism, so the test is computed on the complete
#' cases; impute beforehand if dropping rows is not appropriate for your
#' data. The matrix interface [pe_mediation()] drops nothing: it stops on a
#' missing value and leaves the choice of complete cases to the user.
#'
#' The outcome must be numeric (or logical, for a binary outcome). A factor
#' or character outcome stops with a message, because a factor's numeric
#' codes are not its values: a count stored as a factor would otherwise be
#' analyzed as a different count.
#'
#' @return The tidy `poemed_tbl` returned by [pe_mediation()].
#'
#' @inheritSection POEMED-package How to Cite
#'
#' @author Xiufan Yu and Ken Kelley
#'
#' @seealso [pe_mediation()] for the matrix interface,
#'   [WHO_mediation_design()] for the worked health-expenditure design.
#'
#' @family mediation tests
#'
#' @examples
#' # Build a small data frame and test through the formula interface.
#' set.seed(113)
#' d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
#'                              pattern = "contrasting", c1 = 1)
#' df <- data.frame(y = d$Y, x = d$X[, 1], d$M)
#' meds <- grep("^X", names(df), value = TRUE)   # the mediator columns
#' # The mediator columns are named, so the footer names the selected ones.
#' pe_mediate(y ~ x, data = df, mediators = meds, outcome = "continuous")
#' # HBIC chose the smallest value of the default grid here, so the fit
#' # warned; see the lambda_grid argument of ?pe_mediation.
#'
#' @export
pe_mediate <- function(formula, data,
                       mediators, confounders = NULL,
                       outcome = c("continuous", "binary", "count"), ...) {
  outcome <- .poemed_match_arg(outcome)
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("`formula` must be a two-sided formula, outcome ~ exposure.",
         call. = FALSE)
  if (!is.data.frame(data))
    stop("`data` must be a data.frame.", call. = FALSE)
  # The exposure terms, with a `.` expanded against `data` (every column
  # other than the outcome), so the checks below see the exposures the
  # model will use: y ~ . - m1 - m2 keeps what remains after the
  # subtraction, and a `.` that brings in a mediator or a confounder stops
  # at the role check.
  has_dot <- "." %in% all.vars(formula)
  exposure_terms <- attr(stats::terms(formula, data = data), "term.labels")
  if (!length(exposure_terms))
    stop("`formula` names no exposure; write it as outcome ~ exposure ",
         "(for example y ~ x).", call. = FALSE)

  # The mediators: column names, each listed once, all in `data`, and at
  # least two of them.
  if (!is.character(mediators) || anyNA(mediators))
    stop("`mediators` must be a character vector of mediator column names ",
         "in `data`.", call. = FALSE)
  dup <- unique(mediators[duplicated(mediators)])
  if (length(dup))
    stop("`mediators` lists ", paste(dup, collapse = ", "), " more than ",
         "once; list each mediator column once.", call. = FALSE)
  miss <- setdiff(mediators, names(data))
  if (length(miss))
    stop("mediator column(s) not in `data`: ", paste(miss, collapse = ", "),
         ".", call. = FALSE)
  if (length(mediators) < 2L)
    stop("`mediators` must name at least two mediator columns in `data`; a ",
         "single mediator calls for an ordinary mediation model.",
         call. = FALSE)

  # An empty confounder vector means no confounders, as NULL does.
  if (!length(confounders)) confounders <- NULL
  if (!is.null(confounders)) {
    if (!is.character(confounders) || anyNA(confounders))
      stop("`confounders` must be a character vector of column names in ",
           "`data`, or NULL for none.", call. = FALSE)
    miss <- setdiff(confounders, names(data))
    if (length(miss))
      stop("confounder column(s) not in `data`: ", paste(miss, collapse = ", "),
           ".", call. = FALSE)
  }

  # Each column plays one role. The outcome listed as a mediator would be
  # identified as its own mediator (a p-value of 0), and an exposure or a
  # confounder listed as a mediator changes the test without any sign. The
  # model's columns are the outcome's and those of the exposure terms that
  # remain, so a column only subtracted in `formula` plays no role there.
  outcome_vars <- all.vars(formula[[2L]])
  model_vars <- unique(c(outcome_vars,
                         unlist(lapply(exposure_terms, function(term)
                           all.vars(str2lang(term))))))
  listed <- function(x) {
    paste0(paste(utils::head(x, 10L), collapse = ", "),
           if (length(x) > 10L) sprintf(", and %d more", length(x) - 10L))
  }
  # What to do about a clash: take the column out of the role argument, or,
  # when a `.` brought it in as an exposure, out of the formula.
  fix <- function(clash, role, kind) {
    if (has_dot && length(setdiff(clash, outcome_vars)))
      return(paste0("The `.` in `formula` stands for every column of `data` ",
                    "other than the outcome; name the exposures instead (y ~ ",
                    "x), or subtract the ", kind, " columns (y ~ . - m1 - ",
                    "m2)."))
    paste0("Remove ", if (length(clash) == 1L) "it" else "them", " from `",
           role, "`.")
  }
  clash <- intersect(mediators, model_vars)
  if (length(clash))
    stop("`mediators` must not include the outcome or exposure columns named ",
         "in `formula`: ", listed(clash), ". ",
         fix(clash, "mediators", "mediator"), call. = FALSE)
  clash <- intersect(mediators, confounders)
  if (length(clash))
    stop("`mediators` and `confounders` both list ",
         paste(clash, collapse = ", "), "; list each column in one role only.",
         call. = FALSE)
  clash <- intersect(confounders, model_vars)
  if (length(clash))
    stop("`confounders` must not include the outcome or exposure columns ",
         "named in `formula`: ", listed(clash), ". ",
         fix(clash, "confounders", "confounder"), call. = FALSE)

  # One complete-case rule across every column the model uses. The formula's
  # variables are read with na.pass so that their missing values are counted
  # alongside those of the mediators and confounders; a row missing anything
  # is dropped from all of them, and the user is told how many rows went and
  # where the missing values were. (model.frame()'s own na.omit would drop
  # the formula's incomplete rows silently, and an NA in a factor confounder
  # would then surface as a row-count mismatch in a `Z` the user never built.)
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  used <- c(as.list(mf), as.list(data[mediators]),
            if (!is.null(confounders)) as.list(data[confounders]))
  # A column the formula only subtracts (y ~ . - m1) is still in the model
  # frame; count it once.
  used <- used[!duplicated(names(used))]
  n_missing <- vapply(used, function(v) sum(!stats::complete.cases(v)),
                      integer(1))
  keep <- Reduce(`&`, lapply(used, stats::complete.cases))
  n_drop <- sum(!keep)
  if (n_drop > 0L) {
    cols <- n_missing[n_missing > 0L]
    shown <- utils::head(cols, 10L)
    where <- paste0(paste0(names(shown), " (", shown, ")", collapse = ", "),
                    if (length(cols) > 10L)
                      sprintf(", and %d more columns", length(cols) - 10L)
                    else "")
    if (n_drop == nrow(data))
      stop("Every row of `data` has a missing value in a column the model ",
           "uses (missing values by column: ", where, "); there are no ",
           "complete cases to analyze.", call. = FALSE)
    message(sprintf(paste0(
      "pe_mediate() dropped %d row%s (of %d) with a missing value in a ",
      "column the model uses (missing values by column: %s). POEMED has no ",
      "missing-data mechanism, so the test uses the %d complete cases."),
      n_drop, if (n_drop == 1L) "" else "s", nrow(data), where, sum(keep)))
  }
  mf <- mf[keep, , drop = FALSE]
  data <- data[keep, , drop = FALSE]

  # Outcome and exposures come from the formula (the exposure side may name
  # several columns, e.g. y ~ x1 + x2). A factor or character outcome stops
  # here, where the message can name the column: a factor's numeric codes
  # are not its values.
  Y <- stats::model.response(mf)
  if (is.factor(Y) || is.character(Y)) {
    yname <- paste(deparse(formula[[2L]]), collapse = "")
    stop(sprintf(paste0(
      "The outcome `%s` is a %s, and the test needs its numeric values. ",
      "Convert it first, for example with as.numeric(as.character(%s)) for ",
      "a count or continuous outcome, or code a binary outcome as 0/1."),
      yname, if (is.factor(Y)) "factor" else "character vector", yname),
      call. = FALSE)
  }
  X <- stats::model.matrix(formula, data = mf)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  M <- as.matrix(data[, mediators, drop = FALSE])

  # Confounders: numeric columns pass through, factor/character columns become
  # indicator variables (no intercept column).
  Z <- NULL
  if (!is.null(confounders)) {
    zf <- stats::reformulate(confounders)
    Z <- stats::model.matrix(zf, data = data)
    Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  }

  pe_mediation(X = X, Y = Y, M = M, Z = Z, outcome = outcome, ...)
}
