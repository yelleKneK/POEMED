# Documentation and run-time hygiene that R CMD check reports only with the
# manuals built, or not at all. Ported from the DMAR package and extended to
# every pattern in the portfolio standard's table of what CRAN forbids in
# package code, examples, tests, and vignettes:
#
#   * Rd markup inside \eqn{} or \deqn{} (the PDF manual stops on it);
#   * the two example wrappers the package does without, and an example that
#     runs only when some package is installed;
#   * comment lines in examples that parse as R code;
#   * reading or writing the global environment or the saved RNG state,
#     superassignment, and a negative warn option;
#   * a default file path in a function signature;
#   * more than two cores in examples, tests, or evaluated vignette code, or
#     as a default in package code;
#   * graphics parameters, options, or the working directory changed without
#     a restore on exit;
#   * installing packages, ending the session, or starting external software.
#
# Each detector proves itself on known-bad cases first (a positive control),
# so a detector that cannot fire never passes vacuously. Every name or string
# a detector forbids outright (the example wrappers and guard, the
# global-environment and seed idioms, the warn option, and the names of the
# cores, restore, and forbidden-call checks) is assembled from character
# codes, and so is every known-bad case built on one, so a search of the
# package for any of them finds real offenders only. The markup-in-math and
# commented-code controls are Rd text, which only a help page could carry.

chr <- function(...) intToUtf8(c(...))
rx_escape <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)

tok <- list(
  dontrun   = paste0("\\", chr(100, 111, 110, 116, 114, 117, 110)),
  donttest  = paste0("\\", chr(100, 111, 110, 116, 116, 101, 115, 116)),
  guard     = chr(114, 101, 113, 117, 105, 114, 101, 78, 97, 109, 101, 115, 112, 97, 99, 101),
  rng_state = chr(46, 82, 97, 110, 100, 111, 109, 46, 115, 101, 101, 100),
  global    = chr(46, 71, 108, 111, 98, 97, 108, 69, 110, 118),
  global_fn = chr(103, 108, 111, 98, 97, 108, 101, 110, 118),
  super     = chr(60, 60, 45),
  as_env    = chr(97, 115, 46, 101, 110, 118, 105, 114, 111, 110, 109, 101, 110, 116),
  options   = chr(111, 112, 116, 105, 111, 110, 115),
  warn      = chr(119, 97, 114, 110),
  par       = chr(112, 97, 114),
  setwd     = chr(115, 101, 116, 119, 100),
  detect    = chr(100, 101, 116, 101, 99, 116, 67, 111, 114, 101, 115))

# Calls a package may never make: package installation, ending the session,
# and starting external software (a shell command, a program, a browser).
forbidden_calls <- c(
  chr(105, 110, 115, 116, 97, 108, 108, 46, 112, 97, 99, 107, 97, 103, 101, 115),
  chr(113), chr(113, 117, 105, 116),
  chr(115, 121, 115, 116, 101, 109), chr(115, 121, 115, 116, 101, 109, 50),
  chr(115, 104, 101, 108, 108), chr(115, 104, 101, 108, 108, 46, 101, 120, 101, 99),
  chr(98, 114, 111, 119, 115, 101, 85, 82, 76))

# Argument names that set a number of cores or workers, and the functions
# that start a cluster with the size given first.
core_args <- c("cores", "mc.cores", "ncores", "n_cores", "num_cores", "workers",
               "threads", "nthreads", "n_threads")
cluster_fns <- c("makeCluster", "makePSOCKcluster", "makeForkCluster")

# Rebuilds a function from its source text; the known-bad cases are written
# this way so that their names come from `tok`, never from this file.
fn_from <- function(text) eval(parse(text = text, keep.source = FALSE))
rd_from <- function(text) tools::parse_Rd(textConnection(text))

rd_sources <- function() {
  man <- file.path("..", "..", "man")
  if (dir.exists(man)) tools::Rd_db(dir = file.path("..", "..")) else tools::Rd_db("POEMED")
}

# The Rd walkers recurse and return what they find, so no function here
# assigns outside its own frame.
math_text <- function(x) {
  tag <- attr(x, "Rd_tag")
  here <- if (!is.null(tag) && tag %in% c("\\eqn", "\\deqn")) paste(unlist(x), collapse = "") else character()
  c(here, if (is.list(x)) unlist(lapply(x, math_text)))
}

rd_tags <- function(x) unique(c(attr(x, "Rd_tag"), if (is.list(x)) unlist(lapply(x, rd_tags))))

examples_text <- function(rd) {
  ex <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\examples"), rd)
  if (!length(ex)) return("")
  paste(unlist(ex), collapse = "")
}

markup_in_math <- "\\\\(code|emph|link|strong|pkg|bold|verb|var|samp|file|env|dQuote|sQuote|href|url)\\{"

## ---- Help pages ---------------------------------------------------------------

test_that("the math-markup detector flags the construct CRAN rejected", {
  bad <- rd_from(paste0(
    "\\name{x}\\alias{x}\\title{x}\\description{The region is ",
    "\\eqn{(-\\code{delta_lower}, \\code{delta_upper})}.}"))
  expect_true(any(grepl(markup_in_math, math_text(bad))))
  good <- rd_from(paste0(
    "\\name{x}\\alias{x}\\title{x}\\description{The region is ",
    "\\eqn{(-\\delta_L, +\\delta_U)}, with \\code{delta_lower} as ",
    "\\eqn{\\delta_L}.}"))
  expect_false(any(grepl(markup_in_math, math_text(good))))
})

test_that("no help page carries Rd markup inside \\eqn{} or \\deqn{}", {
  db <- rd_sources()
  skip_if(length(db) == 0L, "no Rd sources available")
  offenders <- names(Filter(function(rd) any(grepl(markup_in_math, math_text(rd))), db))
  expect_identical(offenders, character(0))
})

wrapped_examples <- function(rd) any(c(tok$donttest, tok$dontrun) %in% rd_tags(rd))

test_that("the example-wrapper detector flags both wrappers", {
  for (w in c(tok$donttest, tok$dontrun)) {
    bad <- rd_from(sprintf("\\name{x}\\alias{x}\\title{x}\\description{x}\\examples{\n%s{1 + 1}\n}", w))
    expect_true(wrapped_examples(bad), info = w)
  }
  good <- rd_from("\\name{x}\\alias{x}\\title{x}\\description{x}\\examples{\n1 + 1\n}")
  expect_false(wrapped_examples(good))
})

test_that("no help page wraps its examples", {
  db <- rd_sources()
  skip_if(length(db) == 0L, "no Rd sources available")
  offenders <- names(Filter(wrapped_examples, db))
  expect_identical(offenders, character(0))
})

guarded_examples <- function(rd) grepl(tok$guard, examples_text(rd), fixed = TRUE)

test_that("the example-guard detector flags an example that runs only if a package is present", {
  bad <- rd_from(sprintf(
    "\\name{x}\\alias{x}\\title{x}\\description{x}\\examples{\nif (%s(\"stats\", quietly = TRUE)) 1\n}",
    tok$guard))
  expect_true(guarded_examples(bad))
  good <- rd_from("\\name{x}\\alias{x}\\title{x}\\description{x}\\examples{\nstats::median(1:3)\n}")
  expect_false(guarded_examples(good))
})

test_that("no example block gates on the presence of a package", {
  db <- rd_sources()
  skip_if(length(db) == 0L, "no Rd sources available")
  offenders <- names(Filter(guarded_examples, db))
  expect_identical(offenders, character(0))
})

# Comment lines inside \examples whose text is R code (CRAN's wording:
# "Some code lines in examples are commented out. Please never do that.").
commented_code <- function(rd) {
  ex <- examples_text(rd)
  if (!nzchar(ex)) return(character())
  lines <- strsplit(ex, "\n", fixed = TRUE)[[1]]
  cm <- grep("^\\s*#", lines)
  if (!length(cm)) return(character())
  is_code <- function(txt) {
    parsed <- tryCatch(parse(text = txt, keep.source = FALSE), error = function(e) NULL)
    if (is.null(parsed) || !length(parsed)) return(FALSE)
    is_call <- vapply(parsed, is.call, logical(1))
    is_formula <- vapply(parsed, function(e)
      is.call(e) && identical(e[[1L]], as.name("~")), logical(1))
    any(is_call & !is_formula)
  }
  hits <- integer()
  for (run in split(cm, cumsum(c(1L, diff(cm) != 1L)))) {
    txt <- sub("^\\s*#+\\s?", "", lines[run])
    n <- length(run)
    for (a in seq_len(n)) for (b in a:n) {
      window <- paste(txt[a:b], collapse = "\n")
      if (nzchar(trimws(window)) && is_code(window)) hits <- c(hits, run[a:b])
    }
  }
  lines[sort(unique(hits))]
}

test_that("the commented-code detector separates code from prose", {
  rd <- rd_from(paste0(
    "\\name{x}\\alias{x}\\title{x}\\examples{\n",
    "# The fit is shown for comparison.\n",
    "# pe_mediation(X, Y, M, outcome = \"continuous\")\n",
    "# start <- Sys.time()\n",
    "# Sys.sleep(0.2)\n",
    "mediation_ar1_cov(5, 0.5)\n}"))
  hits <- commented_code(rd)
  expect_length(hits, 3L)
  expect_false(any(grepl("shown for comparison", hits)))
})

test_that("no help page carries commented-out code in its examples", {
  db <- rd_sources()
  skip_if(length(db) == 0L, "no Rd sources available")
  offenders <- Filter(length, lapply(db, commented_code))
  expect_identical(names(offenders), character(0),
                   info = paste(unlist(lapply(names(offenders), function(f)
                     paste(f, offenders[[f]], sep = ": "))), collapse = "\n"))
})

## ---- Package code ---------------------------------------------------------------

namespace_functions <- function() {
  ns <- asNamespace("POEMED")
  nms <- ls(ns, all.names = TRUE)
  fns <- Filter(is.function, mget(nms, envir = ns, inherits = FALSE))
  fns[!vapply(fns, is.primitive, logical(1))]
}

# Every call in an expression (a function's body or its formals), depth
# first. An empty argument, as in x[, 1], is passed over.
all_calls <- function(e) {
  if (!is.call(e) && !is.pairlist(e)) return(list())
  out <- if (is.call(e)) list(e) else list()
  for (i in seq_along(e)) {
    if (identical(e[[i]], quote(expr = ))) next
    out <- c(out, all_calls(e[[i]]))
  }
  out
}
# The name a call invokes, with pkg::fn and pkg:::fn read as fn.
call_name <- function(cl) {
  f <- cl[[1L]]
  if (is.symbol(f)) return(as.character(f))
  if (is.call(f) && length(f) == 3L && as.character(f[[1L]]) %in% c("::", ":::"))
    return(as.character(f[[3L]]))
  ""
}
fn_calls <- function(f) c(all_calls(formals(f)), all_calls(body(f)))
deparse_one <- function(e) paste(deparse(e, width.cutoff = 500L), collapse = " ")

# The global environment, the saved RNG state, superassignment, and a
# negative warn option, read from the deparsed function. Search position 1
# (a pos argument of 1, or the environment at that position) is the global
# environment too.
global_pat <- paste(c(rx_escape(c(tok$global, paste0(tok$global_fn, "()"), tok$rng_state, tok$super)),
                      sprintf("%s\\(%s", tok$options, tok$warn),
                      sprintf("%s\\(1L?\\)", rx_escape(tok$as_env)),
                      sprintf("\\b%s = 1L?\\b", "pos")),
                    collapse = "|")
global_offenders <- function(fns) {
  text <- vapply(fns, function(f) paste(deparse(f, width.cutoff = 500L), collapse = "\n"), character(1))
  names(text)[grepl(global_pat, text)]
}

test_that("the global-environment detector flags every form of a global write or read", {
  bad <- list(
    envir     = fn_from(sprintf("function() assign('z', 1, envir = %s)", tok$global)),
    accessor  = fn_from(sprintf("function() rm(list = 'z', envir = %s())", tok$global_fn)),
    rng       = fn_from(sprintf("function() exists('%s')", tok$rng_state)),
    super     = fn_from(sprintf("function() { k <- 0; f <- function() k %s k + 1; f(); k }", tok$super)),
    warn      = fn_from(sprintf("function() { op <- %s(%s = -1); on.exit(%s(op)); 1 }",
                                tok$options, tok$warn, tok$options)),
    position  = fn_from(sprintf("function() assign('z', 1, %s = 1)", "pos")),
    as_env    = fn_from(sprintf("function() assign('z', 1, envir = %s(1))", tok$as_env)))
  expect_setequal(global_offenders(bad), names(bad))
  good <- list(local = function() { e <- new.env(); assign("z", 1, envir = e); get("z", envir = e) })
  expect_identical(global_offenders(good), character(0))
})

test_that("no function touches the global environment, superassigns, or sets a warn option", {
  expect_identical(global_offenders(namespace_functions()), character(0))
})

# A default file path: any default for a formal named for a file or folder,
# and any default string that ends in a file extension, whatever its formal
# is called.
file_formals <- c("file", "filename", "file_name", "path", "dir", "directory", "outfile")
default_path_offenders <- function(fns) {
  bad <- character()
  for (nm in names(fns)) {
    fm <- formals(fns[[nm]])
    for (arg in names(fm)) {
      if (identical(fm[[arg]], quote(expr = )) || is.null(fm[[arg]])) next
      named_for_file <- arg %in% file_formals
      has_extension <- is.character(fm[[arg]]) && length(fm[[arg]]) == 1L &&
        grepl("^[^[:space:]]*[.][[:alpha:]][[:alnum:]]{0,4}$", fm[[arg]])
      if (named_for_file || has_extension)
        bad <- c(bad, sprintf("%s(%s = %s)", nm, arg, paste(deparse(fm[[arg]]), collapse = "")))
    }
  }
  bad
}

test_that("the default-path detector flags a named formal and any default with a file extension", {
  bad <- list(
    named = fn_from(sprintf("function(x, %s = '%s') NULL", "file", "out.csv")),
    other = fn_from(sprintf("function(x, %s = '%s') NULL", "outfile", "results.csv")))
  expect_length(default_path_offenders(bad), 2L)
  good <- list(ok = function(x, filename = NULL, method = "BH", Y_family = "binomial") NULL)
  expect_identical(default_path_offenders(good), character(0))
})

test_that("no function has a default file path", {
  expect_identical(default_path_offenders(namespace_functions()), character(0))
})

# Graphics parameters, options, and the working directory: a function that
# changes one of them also restores it with on.exit(). A call that only reads
# (par("mfrow"), options("digits")) changes nothing; the working directory
# changes on any call. A change made inside on.exit() is itself the restore.
restore_offenders <- function(fns) {
  bad <- character()
  for (nm in names(fns)) {
    calls <- fn_calls(fns[[nm]])
    names_ <- vapply(calls, call_name, character(1))
    on_exit <- calls[names_ == "on.exit"]
    inside <- unlist(lapply(on_exit, function(oe) vapply(all_calls(oe)[-1L], deparse_one, character(1))))
    restored <- unlist(lapply(on_exit, function(oe) vapply(all_calls(oe)[-1L], call_name, character(1))))
    for (target in c(tok$par, tok$options, tok$setwd)) {
      setters <- Filter(function(cl) {
        if (target == tok$setwd) return(TRUE)
        if (length(cl) < 2L) return(FALSE)
        named <- !is.null(names(cl)) && any(nzchar(names(cl)[-1L]))
        unnamed_value <- any(vapply(seq_along(cl)[-1L], function(i) !is.character(cl[[i]]), logical(1)))
        named || unnamed_value
      }, calls[names_ == target])
      setters <- Filter(function(cl) !deparse_one(cl) %in% inside, setters)
      if (length(setters) && !target %in% restored) bad <- c(bad, sprintf("%s: %s()", nm, target))
    }
  }
  bad
}

test_that("the restore detector flags a change left in place and passes a restored one", {
  bad <- list(
    options = fn_from(sprintf("function() { %s(digits = 3); 1 }", tok$options)),
    par     = fn_from(sprintf("function() { graphics::%s(mfrow = c(1, 2)); 1 }", tok$par)),
    wd      = fn_from(sprintf("function(d) { %s(d); 1 }", tok$setwd)))
  expect_length(restore_offenders(bad), 3L)
  good <- list(
    restored = fn_from(sprintf("function() { op <- %s(mfrow = c(1, 2)); on.exit(%s(op), add = TRUE); 1 }",
                               tok$par, tok$par)),
    reader   = fn_from(sprintf("function() %s(\"digits\")", tok$options)),
    scoped   = function() withr::with_options(list(digits = 3), format(pi)))
  expect_identical(restore_offenders(good), character(0))
})

test_that("no function changes graphics parameters, options, or the working directory without a restore", {
  expect_identical(restore_offenders(namespace_functions()), character(0))
})

forbidden_offenders <- function(fns) {
  bad <- character()
  for (nm in names(fns)) {
    used <- intersect(vapply(fn_calls(fns[[nm]]), call_name, character(1)), forbidden_calls)
    if (length(used)) bad <- c(bad, sprintf("%s: %s()", nm, used))
  }
  bad
}

test_that("the forbidden-call detector flags each forbidden call, namespace-qualified or not", {
  bad <- lapply(forbidden_calls, function(f) fn_from(sprintf("function() { utils::head(1); %s() }", f)))
  names(bad) <- paste0("f", seq_along(bad))
  expect_length(forbidden_offenders(bad), length(forbidden_calls))
  qualified <- list(shell = fn_from(sprintf("function() base::%s()", forbidden_calls[5L])))
  expect_length(forbidden_offenders(qualified), 1L)
  expect_identical(forbidden_offenders(list(ok = function() Sys.time())), character(0))
})

test_that("no function installs packages, ends the session, or starts external software", {
  expect_identical(forbidden_offenders(namespace_functions()), character(0))
})

## ---- Cores --------------------------------------------------------------------

# More than two cores: a core or worker count above 2 given as a number, a
# cluster started with more than two workers, or the machine's core count
# asked for at all.
core_offenders <- function(calls) {
  bad <- character()
  for (cl in calls) {
    nm <- call_name(cl)
    arg_names <- names(cl)
    if (nm == tok$detect) bad <- c(bad, deparse_one(cl))
    for (i in seq_along(cl)[-1L]) {
      big <- is.numeric(cl[[i]]) && length(cl[[i]]) == 1L && cl[[i]] > 2
      named_core <- !is.null(arg_names) && arg_names[i] %in% core_args
      first_cluster_arg <- nm %in% cluster_fns && i == 2L
      if (big && (named_core || first_cluster_arg)) bad <- c(bad, deparse_one(cl))
    }
  }
  unique(bad)
}
code_core_offenders <- function(file, label = basename(file)) {
  exprs <- tryCatch(parse(file, keep.source = FALSE), error = function(e) conditionMessage(e))
  if (is.character(exprs)) return(sprintf("%s: does not parse (%s)", label, exprs))
  hits <- core_offenders(all_calls(as.call(c(as.name("{"), as.list(exprs)))))
  if (length(hits)) sprintf("%s: %s", label, hits) else character()
}
namespace_core_offenders <- function(fns) {
  bad <- character()
  for (nm in names(fns)) {
    fm <- formals(fns[[nm]])
    for (arg in intersect(names(fm), core_args)) {
      if (is.numeric(fm[[arg]]) && length(fm[[arg]]) == 1L && fm[[arg]] > 2)
        bad <- c(bad, sprintf("%s(%s = %s)", nm, arg, fm[[arg]]))
    }
    if (length(core_offenders(fn_calls(fns[[nm]])))) bad <- c(bad, sprintf("%s: %s()", nm, tok$detect))
  }
  unique(bad)
}

test_that("the cores detector flags every way of asking for more than two cores", {
  code <- c(sprintf("pe_power_curve(n = 50, %s = 8L)", "cores"),
            sprintf("%s(%s = 4)", tok$options, "mc.cores"),
            sprintf("parallel::%s(%d)", cluster_fns[1L], 6L),
            sprintf("n <- parallel::%s()", tok$detect),
            sprintf("pe_power_curve(n = 50, %s = 2L)", "cores"))
  file <- tempfile(fileext = ".R")
  on.exit(unlink(file), add = TRUE)
  writeLines(code, file)
  expect_length(code_core_offenders(file), 4L)
  expect_false(any(grepl("= 2L", code_core_offenders(file), fixed = TRUE)))
  bad_default <- list(f = fn_from(sprintf("function(x, %s = 4L) x", "cores")))
  expect_length(namespace_core_offenders(bad_default), 1L)
})

test_that("no function defaults to more than two cores or asks for the core count", {
  expect_identical(namespace_core_offenders(namespace_functions()), character(0))
})

test_that("no example asks for more than two cores", {
  db <- rd_sources()
  skip_if(length(db) == 0L, "no Rd sources available")
  bad <- character()
  for (nm in names(db)) {
    ex <- tempfile(fileext = ".R")
    tools::Rd2ex(db[[nm]], ex)
    if (file.exists(ex)) bad <- c(bad, code_core_offenders(ex, label = nm))
    unlink(ex)
  }
  expect_identical(bad, character(0))
})

test_that("no test asks for more than two cores", {
  files <- list.files(test_path(), pattern = "[.][Rr]$", full.names = TRUE)
  expect_gt(length(files), 1L)
  expect_identical(unlist(lapply(files, code_core_offenders)), character(0))
})

# The evaluated code of every vignette: the sources when they are beside the
# tests (the local suite), purled so that chunks with eval = FALSE come out
# commented; otherwise the code the build purled into the installed doc
# folder (R CMD check). A chunk that is never run may show the article's
# full-scale settings.
test_that("no vignette evaluates code that asks for more than two cores", {
  src <- test_path("..", "..", "vignettes")
  rmd <- if (dir.exists(src)) list.files(src, pattern = "[.]Rmd$", full.names = TRUE) else character()
  if (length(rmd)) {
    skip_if_not_installed("knitr")
    code <- vapply(rmd, function(f) {
      out <- tempfile(fileext = ".R")
      suppressMessages(knitr::purl(f, output = out, quiet = TRUE, documentation = 0L))
    }, character(1))
    names(code) <- basename(rmd)
    on.exit(unlink(code), add = TRUE)
  } else {
    doc <- system.file("doc", package = "POEMED")
    code <- if (nzchar(doc)) list.files(doc, pattern = "[.]R$", full.names = TRUE) else character()
    names(code) <- basename(code)
  }
  skip_if(length(code) == 0L, "no vignette code available")
  bad <- unlist(lapply(names(code), function(v) code_core_offenders(code[[v]], label = v)))
  expect_identical(bad, character(0))
})
