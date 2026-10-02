# Pre-release quality control for POEMED: the static invariants.
#
# Runs every invariant the project conventions state (the package's
# working conventions, and the portfolio QC master where the package does not
# override it), mechanically, and prints one PASS, FAIL, or SKIP line per
# check, with the offending files when a check fails and the reason when a
# check could not run. Exit status is nonzero on any failure, so
# tools/release_gate.R can gate on it; a skipped check is listed, never
# counted as a pass. Modeled on DMAR's tools/release_checks.R; the
# POEMED-specific checks are the citation status of the in-press article,
# the anonymization residue from the double-blind review, the seed idiom,
# the hex logo, the Zenodo metadata, and the pkgdown site configuration.
#
# Usage, from the package root:
#   Rscript tools/release_checks.R
#
# This file is not shipped: tools/ is excluded by .Rbuildignore. Every
# string this file forbids outright (the review-era residue, the run-time
# idioms of the global-environment and seed rules, the em and en dashes and
# a double hyphen used as one, the example blocks the package does without,
# and a truncated package name) is assembled from character codes, and so
# are the known-bad cases the detectors are shown. A search of the
# repository for any of those strings therefore finds real offenders only,
# never this file.

root <- normalizePath(".")
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("Run from the package root.", call. = FALSE)
}

# The tallies live in a small environment, so no function here assigns
# outside its own frame.
tally <- new.env(parent = emptyenv())
tally$failures <- 0L
tally$skipped <- character()
check <- function(name, ok, detail = character()) {
  ok <- isTRUE(ok)
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", name))
  if (!ok) {
    tally$failures <- tally$failures + 1L
    for (d in detail) cat("       ", d, "\n", sep = "")
  }
  invisible(ok)
}
# A check whose input is legitimately absent in this tree is reported as
# skipped, with the reason, and never as passed.
skip_check <- function(name, reason) {
  cat(sprintf("[SKIP] %s\n", name))
  for (d in reason) cat("       ", d, "\n", sep = "")
  tally$skipped <- c(tally$skipped, name)
  invisible(NA)
}

# Character codes to text, and text to a literal regular expression.
chr <- function(...) intToUtf8(c(...))
rx_escape <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)

r_files    <- list.files("R", pattern = "[.]R$", full.names = TRUE)
rd_files   <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
vig_rmd    <- list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE)
prose      <- c(r_files, vig_rmd, "DESCRIPTION", "NEWS.md", "README.md", "inst/CITATION")

grep_files <- function(pattern, files, ...) {
  hits <- character()
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    m <- grep(pattern, lines, ...)
    if (length(m)) hits <- c(hits, sprintf("%s:%d", f, m))
  }
  hits
}
code_lines <- function(f) {
  lines <- readLines(f, warn = FALSE)
  lines[!grepl("^\\s*#", lines)]
}
grep_code <- function(pattern, files = r_files) {
  hits <- character()
  for (f in files) {
    m <- grep(pattern, code_lines(f))
    if (length(m)) hits <- c(hits, sprintf("%s (%d line%s)", f, length(m), if (length(m) == 1) "" else "s"))
  }
  hits
}

## ---- Patterns, and their positive controls -----------------------------
# A detector that cannot fire proves nothing, so each one is first run on a
# known-bad string it must flag (the portfolio's rule 8 for gates). The
# forbidden strings, and the known-bad cases, are built from character codes.
dash2       <- strrep("-", 2L)
dash_pat    <- paste0("[", chr(0x2014, 0x2013), "]")
hyphen_pat  <- sprintf("(^|[^-!<0-9])%s($|[^->0-9a-zA-Z])", dash2)
residue_pat <- paste(c(chr(65, 110, 111, 110, 121, 109, 111, 117, 115),         # the review-era author
                       chr(97, 110, 111, 110, 121, 109, 105, 122, 101, 100),     # "... for review"
                       chr(100, 111, 117, 98, 108, 101, 45, 98, 108, 105, 110, 100),
                       # the pre-acceptance citation states, each in parentheses
                       sprintf("\\(%s\\)", c(chr(109, 105, 110, 111, 114, 32, 114, 101, 118, 105, 115, 105, 111, 110),
                                             chr(115, 117, 98, 109, 105, 116, 116, 101, 100),
                                             chr(117, 110, 100, 101, 114, 32, 114, 101, 118, 105, 101, 119),
                                             chr(105, 110, 32, 112, 114, 101, 112, 97, 114, 97, 116, 105, 111, 110))),
                       chr(99, 111, 110, 100, 105, 116, 105, 111, 110, 97, 108, 108, 121, 32,
                           97, 99, 99, 101, 112, 116, 101, 100),
                       rx_escape(chr(74, 65, 83, 65, 45, 84, 38, 77))),              # the review track
                     collapse = "|")
# The run-time idioms of the seed and global-environment rules.
seed_fn     <- chr(115, 101, 116, 46, 115, 101, 101, 100)                     # the base seed setter
rng_state   <- chr(46, 82, 97, 110, 100, 111, 109, 46, 115, 101, 101, 100)   # the saved RNG state
global_sym  <- chr(46, 71, 108, 111, 98, 97, 108, 69, 110, 118)              # the global environment
global_fn   <- chr(103, 108, 111, 98, 97, 108, 101, 110, 118)                # and its accessor
super_op    <- chr(60, 60, 45)                                                # superassignment
opt_fn      <- chr(111, 112, 116, 105, 111, 110, 115)
warn_arg    <- chr(119, 97, 114, 110)
seed_pat    <- paste0(rx_escape(seed_fn), "\\(")
global_pat  <- paste(rx_escape(c(rng_state, global_sym, paste0(global_fn, "()"))), collapse = "|")
super_pat   <- rx_escape(super_op)
warn_pat    <- sprintf("%s\\(\\s*%s", opt_fn, warn_arg)
# The example blocks the package does without, and the guard it never uses.
dontrun_tag  <- paste0("\\", chr(100, 111, 110, 116, 114, 117, 110))
donttest_tag <- paste0("\\", chr(100, 111, 110, 116, 116, 101, 115, 116))
guard_fn     <- chr(114, 101, 113, 117, 105, 114, 101, 78, 97, 109, 101, 115, 112, 97, 99, 101)
# The package name written short: its first four letters, in any case, with
# no letter or digit before them and not followed by the last two letters
# (so the full name, and every identifier built on it, never matches). An
# underscore before them counts as a separator, so a file name or an
# identifier that carries the short form is caught as well.
short_name  <- chr(80, 79, 69, 77)
short_pat   <- sprintf("(?i)(?<![[:alnum:]])%s(?!ed)", short_name)
is_code <- function(txt) {
  parsed <- tryCatch(parse(text = txt, keep.source = FALSE), error = function(e) NULL)
  if (is.null(parsed) || !length(parsed)) return(FALSE)
  any(vapply(parsed, function(e) is.call(e) && !identical(e[[1L]], as.name("~")), logical(1)))
}
# AP title case: the first and last words and every principal word start
# with a capital; only articles, and conjunctions and prepositions of three
# letters or fewer, stay lowercase. Each part of a hyphenated compound is a
# word. A word written as code (an underscore, a backtick, or a name directly
# followed by an opening parenthesis, as in a call or AR(1)) is left as is; a
# parenthesis that only opens or closes a parenthetical is punctuation, so the
# words inside a parenthetical are checked like any other.
ap_small <- c("a", "an", "the", "and", "but", "for", "nor", "or", "so", "yet",
              "as", "at", "by", "in", "of", "off", "on", "out", "per", "to", "up", "via", "vs")
title_case_problems <- function(titles) {
  bad <- vapply(titles, function(t) {
    words <- strsplit(trimws(t), "\\s+")[[1]]
    words <- words[!grepl("[_`]|[[:alnum:]]\\(", words)]
    parts <- unlist(strsplit(words, "-", fixed = TRUE))
    parts <- gsub("^[^[:alpha:]]+|[^[:alpha:]]+$", "", parts)
    parts <- parts[nzchar(parts)]
    if (!length(parts)) return(FALSE)
    lower <- grepl("^[[:lower:]]", parts)
    edge <- c(1L, length(parts))
    any(lower[edge]) || any(lower[-edge] & !tolower(parts[-edge]) %in% ap_small)
  }, logical(1))
  titles[bad]
}
# The headings of a help page: its title and the title of each of its
# sections, as text. Code-like markup (\code{}, \pkg{}, and the like) is
# written between backticks, with any space inside it replaced, so that a
# package or function name keeps its own case and is not read as a word.
rd_code_tags <- paste0("\\", c("code", "pkg", "link", "var", "env", "file", "samp",
                               "verb", "eqn", "deqn", "href", "url"))
rd_text <- function(x) {
  tag <- attr(x, "Rd_tag")
  if (!is.null(tag) && tag %in% rd_code_tags)
    return(paste0("`", gsub("\\s+", "_", paste(unlist(x), collapse = "")), "`"))
  if (is.list(x)) return(paste(vapply(x, rd_text, ""), collapse = ""))
  paste(as.character(x), collapse = "")
}
rd_headings <- function(rd) {
  tags <- vapply(rd, function(x) if (is.null(attr(x, "Rd_tag"))) "" else attr(x, "Rd_tag"), "")
  heads <- c(title = unlist(lapply(rd[tags == "\\title"], rd_text)),
             section = unlist(lapply(rd[tags == "\\section"], function(x) rd_text(x[[1L]]))))
  gsub("\\s+", " ", trimws(heads))
}
# The formulas of a help page: how many \eqn{} and \deqn{} it has, and the
# LaTeX of each one with a single argument, which the text help (?topic in a
# terminal) shows raw. The second argument is the plain-text form.
rd_math <- function(x) {
  tag <- attr(x, "Rd_tag")
  is_math <- !is.null(tag) && tag %in% c("\\eqn", "\\deqn")
  inner <- if (is.list(x) && !is_math) lapply(x, rd_math) else list()
  list(n = is_math + sum(vapply(inner, `[[`, 0, "n")),
       bare = c(if (is_math && length(x) < 2L)
                  gsub("\\s+", " ", trimws(paste(unlist(x), collapse = ""))),
                unlist(lapply(inner, `[[`, "bare"))))
}
# The Articles menu of the pkgdown site: every article other than the
# package's own vignette (the "Get started" link) sits in a section that
# carries a navbar key, or pkgdown leaves it out of the menu.
articles_off_menu <- function(cfg, vignettes, package) {
  on_menu <- unlist(lapply(cfg$articles, function(s) if ("navbar" %in% names(s)) unlist(s$contents)))
  setdiff(setdiff(vignettes, package), on_menu)
}
controls <- c(
  dash      = grepl(dash_pat, paste0("a ", chr(0x2014), " b")),
  hyphen    = grepl(hyphen_pat, paste("the benchmark", dash2, "the one")) &&
    !grepl(hyphen_pat, paste0("2000", dash2, "2021")),
  residue   = grepl(residue_pat, sprintf("Yu and Kelley (%s)", chr(109, 105, 110, 111, 114, 32, 114, 101, 118, 105, 115, 105, 111, 110)),
                    ignore.case = TRUE),
  seed      = grepl(seed_pat, paste0(seed_fn, "(1)")) && !grepl(seed_pat, "withr::local_seed(1)"),
  global    = all(grepl(global_pat, c(sprintf("exists('%s')", rng_state), paste0("assign('z', 1, envir = ", global_sym, ")"),
                                      paste0(global_fn, "()")))),
  superassign = grepl(super_pat, paste("k", super_op, "k + 1")) && !grepl(super_pat, "k <- k + 1"),
  warn      = grepl(warn_pat, sprintf("%s(%s = -1)", opt_fn, warn_arg)),
  examples  = grepl(rx_escape(donttest_tag), paste0(donttest_tag, "{x}"), perl = TRUE) &&
    grepl(rx_escape(dontrun_tag), paste0(dontrun_tag, "{x}"), perl = TRUE),
  code      = is_code("fit <- pe_mediation(X, Y, M)") && !is_code("the mediators all push the same way"),
  short_name = all(grepl(short_pat, c(paste("the", short_name, "tests"), paste0("hex_", short_name, ".svg"),
                                      paste0(tolower(short_name), "_1.0.0.tar.gz")), perl = TRUE)) &&
    !any(grepl(short_pat, c("POEMED", "poemed_tbl", "hex_POEMED.svg", "POEMED's"), perl = TRUE)),
  title_case = identical(unname(title_case_problems(c("Testing for mediation", "Testing for Mediation",
                                                      "Results, Printing, and Tidiers", "A Head-to-Head",
                                                      "The Article and the Benchmark",
                                                      "Mediation Data (benchmark data set)",
                                                      "Autoregressive (AR(1)) Covariance Matrix",
                                                      "Working With `broom`"))),
                         c("Testing for mediation", "Mediation Data (benchmark data set)")),
  rd_headings = local({
    rd <- tools::parse_Rd(textConnection(c(
      "\\name{x}", "\\alias{x}", "\\title{Broom verbs for \\pkg{POEMED} results}",
      "\\description{d}", "\\section{Subsetting and combining}{s}",
      "\\section{Working With \\pkg{broom}}{s}")))
    identical(unname(title_case_problems(rd_headings(rd))),
              c("Broom verbs for `POEMED` results", "Subsetting and combining"))
  }),
  rd_math = identical(rd_math(tools::parse_Rd(textConnection(c(
    "\\name{x}", "\\alias{x}", "\\title{X}",
    "\\description{\\eqn{\\beta} and \\eqn{\\alpha}{alpha}, \\deqn{J_m}{J_m}}")))),
    list(n = 3, bare = "\\beta")),
  articles_menu = identical(articles_off_menu(list(articles = list(
    list(title = "Start Here", navbar = NULL, contents = list("pkg")),
    list(title = "More", contents = list("one", "two")))), c("pkg", "one", "two"), "pkg"), c("one", "two")))
check("detector positive controls (each detector flags its known-bad case)",
      all(controls), names(controls)[!controls])

## ---- Examples ---------------------------------------------------------------

h <- grep_files(rx_escape(dontrun_tag), c(r_files, rd_files))
check(sprintf("no %s anywhere", dontrun_tag), length(h) == 0, h)

# No example block of the other kind anywhere (the portfolio master, decided
# 2026-08-10 and adopted for POEMED on 2026-09-25): one such block makes the
# CRAN-mode check run the whole example corpus twice, so every example is
# kept fast enough to run.
h <- grep_files(rx_escape(donttest_tag), c(r_files, rd_files))
check(sprintf("no %s anywhere", donttest_tag), length(h) == 0, h)

# No commented-out code in any @examples block (CRAN: "Some code lines in
# examples are commented out. Please never do that."). Each contiguous
# window of comment lines is parsed with the markers stripped; prose does
# not parse and is never flagged.
commented_hits <- character()
for (f in r_files) {
  lines <- readLines(f, warn = FALSE)
  for (s in grep("^#' @examples", lines)) {
    i <- s + 1L; ex <- character(); idx <- integer()
    while (i <= length(lines) && grepl("^#'", lines[i]) && !grepl("^#' @[a-zA-Z]", lines[i])) {
      ex <- c(ex, sub("^#' ?", "", lines[i])); idx <- c(idx, i); i <- i + 1L
    }
    cm <- grep("^\\s*#", ex)
    if (!length(cm)) next
    for (run in split(cm, cumsum(c(1L, diff(cm) != 1L)))) {
      txt <- sub("^\\s*#+\\s?", "", ex[run]); n <- length(run); flagged <- integer()
      for (a in seq_len(n)) for (b in a:n) {
        w <- paste(txt[a:b], collapse = "\n")
        if (nzchar(trimws(w)) && is_code(w)) flagged <- c(flagged, run[a:b])
      }
      if (length(flagged)) commented_hits <- c(commented_hits, sprintf("%s:%d", f, idx[sort(unique(flagged))]))
    }
  }
}
check("no commented-out code in any @examples block", length(commented_hits) == 0, commented_hits)

ex_lines <- function(f) {
  lines <- readLines(f, warn = FALSE); out <- integer(); in_ex <- FALSE
  for (i in seq_along(lines)) {
    if (grepl("^#' @examples", lines[i])) { in_ex <- TRUE; next }
    if (in_ex && (!grepl("^#'", lines[i]) || grepl("^#' @[a-zA-Z]", lines[i]))) in_ex <- FALSE
    if (in_ex) out <- c(out, i)
  }
  out
}
guard_hits <- character(); seed_hits <- character()
for (f in r_files) {
  lines <- readLines(f, warn = FALSE)
  for (i in ex_lines(f)) {
    if (grepl(guard_fn, lines[i], fixed = TRUE) && grepl("^#' [^#]", lines[i]))
      guard_hits <- c(guard_hits, sprintf("%s:%d", f, i))
    if (grepl(seed_pat, lines[i]) && !grepl(paste0(seed_pat, "113\\)"), lines[i]))
      seed_hits <- c(seed_hits, sprintf("%s:%d", f, i))
  }
}
check(sprintf("no %s() guard executing in @examples", guard_fn), length(guard_hits) == 0, guard_hits)
check(sprintf("every example seed is %s(113)", seed_fn), length(seed_hits) == 0, seed_hits)

## ---- Run-time policies (code lines only) -------------------------------------

# Any function using randomness takes seed = NULL and goes through
# .poemed_local_seed() (withr::local_seed underneath), which sets the seed for
# the call and restores the caller's state on exit. Nothing in R/ calls the
# base seed setter or reads or writes the saved RNG state or the global
# environment: that idiom is what CRAN's review of DMAR named on 2026-09-04.
h <- grep_code(seed_pat)
check(sprintf("no %s() anywhere in R/ (seeds go through .poemed_local_seed())", seed_fn), length(h) == 0, h)
h <- grep_code(global_pat)
check(sprintf("no function reads or writes %s, %s, or %s()", rng_state, global_sym, global_fn), length(h) == 0, h)
seed_default <- grep_files("^[^#]*function\\(.*seed *= *[0-9]", r_files)
check("no baked-in default seed in any signature", length(seed_default) == 0, seed_default)
h <- grep_code(super_pat)
check(sprintf("no %s anywhere in R/", super_op), length(h) == 0, h)
h <- grep_code(warn_pat)
check(sprintf("no %s(%s = ...) anywhere in R/", opt_fn, warn_arg), length(h) == 0, h)
h <- grep_code("(filename|file|path)\\s*=\\s*[\"']")
check("no default file path in any function signature or call", length(h) == 0, h)

# Full precision: round(), signif(), and formatC() only in the display layer
# or in an integer-validation comparison.
display <- c("R/poemed_tbl.R", "R/summary_poemed_tbl.R")
rnd <- character()
for (f in setdiff(r_files, display)) {
  cl <- code_lines(f)
  m <- grep("\\b(round|signif|formatC)\\(", cl)
  m <- m[!grepl("(==|!=|-)\\s*round\\(", cl[m])]
  if (length(m)) rnd <- c(rnd, sprintf("%s: %s", f, trimws(cl[m])))
}
check("no round()/signif()/formatC() outside the display layer and validation",
      length(rnd) == 0, rnd)

## ---- DESCRIPTION, NEWS, citation -----------------------------------------------

desc <- read.dcf("DESCRIPTION")
ver <- unname(desc[, "Version"])
pkg_name <- unname(desc[, "Package"])
squish <- function(x) gsub("\\s+", " ", trimws(unname(x)))
news_head <- grep("^# ", readLines("NEWS.md", warn = FALSE), value = TRUE)[1]
check(sprintf("NEWS.md leads with version %s", ver),
      identical(sub("^# [A-Za-z]+ ", "", news_head), ver), news_head)

auth <- eval(parse(text = desc[, "Authors@R"]))
roles <- lapply(auth, function(p) p$role)
fam <- vapply(auth, function(p) p$family, "")
auth_ok <- identical(fam, c("Yu", "Kelley")) &&
  identical(roles[[1]], "aut") && setequal(roles[[2]], c("aut", "cre")) &&
  identical(unname(auth[[2]]$email), "kkelley@nd.edu")
check("Authors@R: Yu (aut), Kelley (aut, cre, kkelley@nd.edu)", auth_ok,
      format(auth, include = c("given", "family", "role", "email")))
url_ok <- "URL" %in% colnames(desc) && grepl("github.com/yelleKneK/POEMED", desc[, "URL"]) &&
  "BugReports" %in% colnames(desc) && grepl("github.com/yelleKneK/POEMED/issues", desc[, "BugReports"])
check("URL and BugReports name github.com/yelleKneK/POEMED", url_ok)

# The article's citation status is one fact stated in several places, and
# the places must agree: while inst/CITATION carries note = "In press", the
# DESCRIPTION and README cite it as in press; once it carries a volume, no
# "in press" may remain anywhere.
cit_txt <- paste(readLines("inst/CITATION", warn = FALSE), collapse = "\n")
published <- grepl("\\bvolume\\s*=", cit_txt)
in_press_hits <- grep_files("in press", c(prose, man = rd_files), ignore.case = TRUE)
if (published) {
  check("article published: no 'in press' left anywhere", length(in_press_hits) == 0, in_press_hits)
} else {
  need <- c(inst_CITATION = grepl('note\\s*=\\s*"In press"', cit_txt),
            DESCRIPTION = grepl("(in press)", desc[, "Description"], fixed = TRUE),
            README = any(grepl("(in press)", readLines("README.md", warn = FALSE), fixed = TRUE)),
            CITATION.cff = !file.exists("CITATION.cff") ||
              any(grepl("status: in-press", readLines("CITATION.cff", warn = FALSE), fixed = TRUE)))
  check("article in press: CITATION, DESCRIPTION, README, and CITATION.cff agree",
        all(need), names(need)[!need])
}
h <- grep_files(residue_pat, c(prose, rd_files), ignore.case = TRUE)
check("no review-era residue (anonymized authors, pre-acceptance citation status)",
      length(h) == 0, h)

cit <- tryCatch(utils::readCitationFile("inst/CITATION", meta = as.list(desc[1, ])),
                error = function(e) conditionMessage(e))
check("inst/CITATION readable from DESCRIPTION metadata alone (CRAN incoming state)",
      inherits(cit, "bibentry") && !grepl("packageVersion", cit_txt),
      if (!inherits(cit, "bibentry")) cit)

lock <- c(CITATION.cff = !file.exists("CITATION.cff") ||
            any(grepl(sprintf("^version: %s$", ver), readLines("CITATION.cff", warn = FALSE))))
check(sprintf("version %s in lockstep (DESCRIPTION, NEWS, CITATION.cff)", ver),
      all(lock), names(lock)[!lock])
# The working-conventions file (untracked; its name is assembled from
# character codes so this tracked script never carries it) states the
# version in its status block and must agree. A tree without it (a fresh
# clone) reports the check as skipped, not passed.
conv_file <- chr(67, 76, 65, 85, 68, 69, 46, 109, 100)
if (file.exists(conv_file)) {
  check(sprintf("version %s on the conventions status line", ver),
        any(grepl(sprintf("Version: %s", ver), readLines(conv_file, warn = FALSE), fixed = TRUE)),
        sprintf("the conventions file has no line reading Version: %s", ver))
} else {
  skip_check(sprintf("version %s on the conventions status line", ver),
             "the untracked conventions file is not in this tree")
}

# Zenodo builds its record from .zenodo.json at each GitHub release, so the
# record says exactly what DESCRIPTION says: the same title (after the
# package name), the same version, and DESCRIPTION's Description as its one
# HTML paragraph, compared after whitespace is normalized.
zen <- tryCatch(jsonlite::read_json(".zenodo.json"), error = function(e) conditionMessage(e))
if (is.character(zen)) {
  check(".zenodo.json reads as JSON", FALSE, zen)
} else {
  zen_ok <- c(title = identical(zen$title, sprintf("%s: %s", pkg_name, squish(desc[, "Title"]))),
              version = identical(zen$version, ver),
              description = identical(squish(gsub("</?p>", "", zen$description)),
                                      squish(desc[, "Description"])))
  check(".zenodo.json title, version, and description agree with DESCRIPTION",
        all(zen_ok), names(zen_ok)[!zen_ok])
}

## ---- House style ---------------------------------------------------------------

h <- grep_files(dash_pat, prose)
check("no em or en dashes", length(h) == 0, h)
h <- grep_files(hyphen_pat, setdiff(prose, "DESCRIPTION"))
h <- h[!grepl("R CMD|--as-cran|--run-donttest", vapply(h, function(x) {
  p <- strsplit(x, ":")[[1]]; readLines(p[1], warn = FALSE)[as.integer(p[2])] }, ""))]
check("no double hyphen used as a dash in prose (numeric ranges excepted)", length(h) == 0, h)
nonascii <- character()
for (f in c(r_files, "DESCRIPTION")) {
  m <- grep("[^\\x01-\\x7F]", readLines(f, warn = FALSE, encoding = "UTF-8"), perl = TRUE)
  if (length(m)) nonascii <- c(nonascii, sprintf("%s:%d", f, m))
}
check("no non-ASCII characters in R/ or DESCRIPTION", length(nonascii) == 0, nonascii)

## ---- Structure, data, vignettes, logo ----------------------------------------

data_objs <- sub("[.]rda$", "", list.files("data", pattern = "[.]rda$"))
rd_db <- tools::Rd_db(dir = ".")
aliases <- unlist(lapply(rd_db, function(rd) {
  unlist(lapply(rd[vapply(rd, function(x) identical(attr(x, "Rd_tag"), "\\alias"), logical(1))],
                function(x) paste(unlist(x), collapse = "")))
}))
has_fs <- vapply(rd_db, function(rd) {
  tags <- vapply(rd, function(x) attr(x, "Rd_tag"), "")
  all(c("\\format", "\\source") %in% tags)
}, logical(1))
fs_aliases <- unlist(lapply(names(rd_db)[has_fs], function(n) {
  rd <- rd_db[[n]]
  unlist(lapply(rd[vapply(rd, function(x) identical(attr(x, "Rd_tag"), "\\alias"), logical(1))],
                function(x) paste(unlist(x), collapse = "")))
}))
# example_* objects share one page (example_data.Rd).
undoc <- data_objs[!data_objs %in% fs_aliases]
check("every shipped data set documented with \\format and \\source", length(undoc) == 0, undoc)

# Every help-page title and section heading is in AP title case, like the
# vignette and pkgdown headings (the conventions cover headings and titles).
rd_heads <- unlist(lapply(names(rd_db), function(n) {
  h <- rd_headings(rd_db[[n]])
  stats::setNames(h, rep(n, length(h)))
}))
bad_heads <- title_case_problems(rd_heads)
check(sprintf("help-page titles and section headings in AP title case (%d headings)", length(rd_heads)),
      length(rd_heads) > 0 && !length(bad_heads), sprintf("%s: %s", names(bad_heads), bad_heads))

# Every formula on a help page carries a plain-text form as its second
# argument, so the text help shows no raw LaTeX.
math <- lapply(rd_db, rd_math)
no_ascii <- unlist(lapply(names(math), function(n)
  if (length(math[[n]]$bare)) sprintf("%s: %s", n, math[[n]]$bare)))
n_math <- sum(vapply(math, `[[`, 0, "n"))
check(sprintf("every \\eqn{} and \\deqn{} has a plain-text second argument (%d formulas)", n_math),
      n_math > 0 && !length(no_ascii), no_ascii)

vig_names <- sub("[.]Rmd$", "", basename(vig_rmd))
dangling <- character()
for (hit in grep_files('vignette\\("([^"]+)"', c(r_files, vig_rmd, "README.md"))) {
  parts <- strsplit(hit, ":")[[1]]
  line <- readLines(parts[1], warn = FALSE)[as.integer(parts[2])]
  for (call in regmatches(line, gregexpr('vignette\\("([^"]+)"', line))[[1]]) {
    nm <- sub('vignette\\("([^"]+)".*', "\\1", call)
    if (!nm %in% vig_names) dangling <- c(dangling, paste0(hit, " -> ", nm))
  }
}
check("no vignette() call names a vignette that does not ship", length(dangling) == 0, dangling)

# The pkgdown site. pkgdown writes the help pages' \eqn{} and \deqn{} as TeX
# and renders it only when the template declares KaTeX or MathJax (its
# default, MathML, reaches the articles alone). The Articles menu lists every
# article other than the package's own vignette, which is the "Get started"
# link. Every section title is in AP title case.
cfg <- tryCatch(yaml::read_yaml("_pkgdown.yml"), error = function(e) conditionMessage(e))
if (is.character(cfg)) {
  check("_pkgdown.yml reads as YAML", FALSE, cfg)
} else {
  math_pages <- basename(rd_files)[vapply(rd_files, function(f)
    any(grepl("\\\\(eqn|deqn)\\{", readLines(f, warn = FALSE))), logical(1))]
  renderer <- cfg$template$`math-rendering`
  check(sprintf("pkgdown renders the formulas on %d help pages (template math-rendering is katex or mathjax)",
                length(math_pages)),
        !length(math_pages) || isTRUE(renderer %in% c("katex", "mathjax")),
        sprintf("template math-rendering is %s", if (is.null(renderer)) "unset (MathML, which leaves raw TeX)" else renderer))
  off_menu <- articles_off_menu(cfg, vig_names, pkg_name)
  check("every article other than the package's own vignette is in the pkgdown Articles menu",
        !length(off_menu), sprintf("not in the menu: %s", off_menu))
  sec_titles <- unlist(lapply(c(cfg$reference, cfg$articles), function(sec) sec$title))
  bad_titles <- title_case_problems(sec_titles)
  check(sprintf("pkgdown section titles in AP title case (%d titles)", length(sec_titles)),
        length(sec_titles) > 0 && !length(bad_titles), bad_titles)
}

# A git that actually runs. On this machine /usr/local/bin/git is a stale
# Intel binary that precedes /usr/bin in the PATH R's shell sees; it exits
# 126, and a git check fed its empty output would pass having checked
# nothing. So resolve a working git first and fail loudly without one.
git_bin <- local({
  cands <- unique(c(Sys.which("git"), "/usr/bin/git", "/opt/homebrew/bin/git"))
  cands <- cands[nzchar(cands) & file.exists(cands)]
  ok <- vapply(cands, function(g) identical(suppressWarnings(
    tryCatch(system2(g, "--version", stdout = FALSE, stderr = FALSE), error = function(e) 1L)), 0L), logical(1))
  if (any(ok)) cands[ok][1] else NA_character_
})
tracked <- if (is.na(git_bin)) character() else
  suppressWarnings(system2(git_bin, "ls-files", stdout = TRUE))
check("a working git is available (the git checks below are not vacuous)",
      !is.na(git_bin) && length(tracked) > 0 && is.null(attr(tracked, "status")),
      "no git binary on this machine ran; checked Sys.which, /usr/bin, /opt/homebrew/bin")
gen <- grep("^vignettes/.*[.](R|html)$|^doc/|^Meta/|^inst/doc/|[.]Rcheck/|[.]tar[.]gz$|Rplots[.]pdf$",
            tracked, value = TRUE)
check("no generated build products tracked in git", length(gen) == 0, gen)
private <- grep(paste0("(^|/)", conv_file, "$|[.]docx$|GAP-ANALYSIS|HANDOFF|QC_|_QC|RELEASE_QC|evaluation|audit"),
                tracked, value = TRUE, ignore.case = TRUE)
check("no private working documents tracked in git", length(private) == 0, private)

# The package name is written in full everywhere: no tracked text file
# carries its first four letters standing alone (short_pat, above). Binary
# files (a NUL byte in the first 8000 bytes, git's own test) are passed over,
# and data URIs are stripped first so that encoded bytes cannot match.
is_text_file <- function(f) {
  if (!file.exists(f) || dir.exists(f)) return(FALSE)
  con <- file(f, "rb")
  on.exit(close(con))
  !any(readBin(con, "raw", 8000L) == as.raw(0L))
}
text_files <- tracked[vapply(tracked, is_text_file, logical(1))]
short_hits <- character()
for (f in text_files) {
  lines <- gsub("data:[^\"')[:space:]]+", "", suppressWarnings(readLines(f, warn = FALSE)), useBytes = TRUE)
  m <- grep(short_pat, lines, perl = TRUE, useBytes = TRUE)
  if (length(m)) short_hits <- c(short_hits, sprintf("%s:%d", f, m))
}
check(sprintf("the package name is written in full in all %d tracked text files (no truncated form)",
              length(text_files)),
      length(text_files) > 0 && !length(short_hits), short_hits)

# The hex logo: shipped, referenced, and still the family master's art. The
# family master and the generator are named from DESCRIPTION's Package
# field, and a master or generator that is missing fails its check: a
# comparison that cannot run is not a pass.
logo <- c(png = file.exists("man/figures/logo.png"), svg = file.exists("man/figures/logo.svg"),
          readme = any(grepl("man/figures/logo.png", readLines("README.md", warn = FALSE), fixed = TRUE)))
check("hex logo shipped (man/figures/logo.png, logo.svg) and shown in README", all(logo),
      names(logo)[!logo])
same_bytes <- function(a, b) file.exists(a) && file.exists(b) &&
  identical(unname(tools::md5sum(a)), unname(tools::md5sum(b)))
family_svg <- file.path(dirname(root), "KenKelleyHexFamily", "hexes", sprintf("hex_%s.svg", pkg_name))
check(sprintf("logo.svg identical to the KenKelleyHexFamily master (hexes/%s)", basename(family_svg)),
      same_bytes("man/figures/logo.svg", family_svg),
      c(if (!file.exists(family_svg)) sprintf("family master not found: %s", family_svg),
        if (!logo[["svg"]]) "man/figures/logo.svg not found"))
hex_gen <- file.path("dev", sprintf("make_%s_hex.R", pkg_name))
regen_ok <- FALSE
regen_detail <- character()
if (!file.exists(hex_gen)) {
  regen_detail <- sprintf("generator not found: %s", hex_gen)
} else {
  gen_dir <- tempfile("hex")
  dir.create(file.path(gen_dir, "dev"), recursive = TRUE)
  file.copy(hex_gen, file.path(gen_dir, "dev"))
  st <- withr::with_dir(gen_dir, system2("Rscript", hex_gen, stdout = FALSE, stderr = FALSE))
  regen_ok <- identical(as.integer(st), 0L) &&
    same_bytes(file.path(gen_dir, "man", "figures", "logo.svg"), "man/figures/logo.svg")
  if (!identical(as.integer(st), 0L)) regen_detail <- sprintf("%s exited with status %s", hex_gen, st)
}
check(sprintf("%s regenerates the shipped logo.svg byte for byte", hex_gen), regen_ok, regen_detail)

# Tarball hygiene: build one (vignettes off, so the size is a floor) and
# confirm nothing ships that should not.
cat("       building the tarball to inspect it...\n")
bld <- tempfile("bld"); dir.create(bld)
old_wd <- setwd(bld)
ok_build <- system2("R", c("CMD", "build", "--no-manual", "--no-build-vignettes", shQuote(root)),
                    stdout = FALSE, stderr = FALSE) == 0
setwd(old_wd)
if (ok_build) {
  tb <- list.files(bld, pattern = "[.]tar[.]gz$", full.names = TRUE)[1]
  contents <- untar(tb, list = TRUE)
  bad <- c(grep("tools/|dev/|[.]github|_pkgdown|CITATION[.]cff|[.]Rcheck|cran-comments|logo-hires|Rplots",
                contents, value = TRUE),
           grep("^[^/]+/(?!NEWS[.]md$|README[.]md$|LICENSE[.]md$)[^/]+[.]md$",
                contents, value = TRUE, perl = TRUE))
  check("tarball contains no working files", length(bad) == 0, bad)
  sz <- file.info(tb)$size / 1048576
  check(sprintf("tarball floor %.2f Mb below the 5 Mb guideline", sz), sz < 5)
} else {
  check("R CMD build succeeds", FALSE, "build failed; run it by hand for the error")
}

cat(sprintf("\n%d failure(s), %d skipped%s.\n", tally$failures, length(tally$skipped),
            if (length(tally$skipped)) paste0(" (", paste(tally$skipped, collapse = "; "), ")") else ""))
if (tally$failures > 0) quit(status = 1L)
