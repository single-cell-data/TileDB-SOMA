
#' Signal a Deprecation
#'
#' Signal a deprecation/defunct stage using \CRANpkg{lifecycle}
#'
#' @inherit .deprecation_stage details
#'
#' @inheritParams lifecycle::deprecate_warn
#' @param ... Must be empty, and error will be thrown if anything is passed
#' through the \dots
#' @param env,user_env Needed for \
#' code{\link[lifecycle:deprecate_warn]{lifecycle::deprecate_*}()},
#' do not modify
#'
#' @return No return value, used for the side effects of signaling a
#' deprecation or defunct stage
#'
#' @keywords internal
#'
#' @seealso \code{\link[lifecycle:deprecate_warn]{lifecycle::deprecate_warn}()},
#' \code{\link[lifecycle:deprecate_stop]{lifecycle::deprecate_stop}()},
#' \code{\link{.deprecation_stage}()}
#'
#' @noRd
#'
#' @examples
#' \dontrun{\donttest{
#' SOMANDArrayBase <- R6::R6Class(
#'   "SOMANDArrayBase",
#'   ...,
#'   public = list(
#'     set_data_type(type) {
#'       .deprecate(when = "2.0.0", what = "SOMANDArrayBase$set_data_type()")
#'     }
#'   )
#' )
#' }}
#'
.deprecate <- function(
  when,
  what,
  with = NULL,
  ...,
  details = NULL,
  id = NULL,
  always = TRUE,
  env = rlang::caller_env(),
  user_env = rlang::caller_env(n = 2L)
) {
  stopifnot(
    rlang::is_character(what, n = 1L),
    is.null(x = with) || rlang::is_character(with)
  )
  rlang::check_dots_empty0(...)
  switch(
    EXPR = .deprecation_stage(when = when) %||% "future",
    defunct = lifecycle::deprecate_stop(
      when = as.character(when),
      what = I(what),
      with = with,
      details = details,
      env = env
    ),
    deprecate = lifecycle::deprecate_warn(
      when = as.character(when),
      what = I(what),
      with = with,
      details = details,
      # lifecycle tries to be clever when determining when to warn; however,
      # it's actually pretty bad at it. It doesn't work well with R6 nor does
      # it accurately track testthat usage. We need to force it
      # to throw a deprecation warning
      id = id %||% as.character(Sys.time()),
      always = always,
      env = env,
      user_env = user_env
    )
  )
}

#' Determine the Deprecation Stage
#'
#' Attempts to determine the appropriate deprecation stage following
#' TileDB-SOMA's
#' \href{https://github.com/single-cell-data/TileDB-SOMA/blob/main/dev_docs/POLICIES.md}{deprecation policy}.
#'
#' \pkg{tiledbsoma} follows three stages of deprecation. The three stages are
#' \itemize{
#'  \item no deprecation
#'  \item deprecated, signaled with a \link[lifecycle:deprecate_warn]{warning}
#'  \item defunct, signaled with an \link[lifecycle:deprecate_stop]{error}
#' }
#' The defunct stage is signaled if:
#' \itemize{
#'  \item the major version of \pkg{tiledbsoma} is greater than \code{when}, or
#'  \item the major versions of \pkg{tiledbsoma} and \code{when} are the same,
#'   the minor version of \pkg{tiledbsoma} is higher than \code{when} by two or
#'   more, and it has been at least three months (twelve weeks) since
#'   \code{when} was released
#' }
#' The deprecated stage is signaled if the major versions of \pkg{tiledbsoma}
#' and \code{when} are the same and
#' \itemize{
#'  \item the minor version of \pkg{tiledbsoma} is the same or only one higher
#'   than \code{when}, \strong{or}
#'  \item it has been less than three months (twelve weeks) since \code{when}
#'   was released
#' }
#' If the major version of \pkg{tiledbsoma} is less than \code{when}, or if the
#' major versions are the same, but minor version of \pkg{tiledbsoma} is less
#' than \code{when}, no deprecation is signaled
#'
#' @inheritParams lifecycle::deprecate_warn
#'
#' @return Depending on the appropriate deprecation stage, returns one of:
#' \itemize{
#'  \item \dQuote{\code{deprecate}} for the deprecated stage
#'  \item \dQuote{\code{defunct}} for the defunct stage
#'  \item \code{NULL} invisibly for the no deprecation stage
#' }
#'
#' @keywords internal
#'
#' @seealso \code{\link[lifecycle:deprecate_warn]{lifecycle::deprecate_warn}()},
#' \code{\link[lifecycle:deprecate_stop]{lifecycle::deprecate_stop}()},
#' \code{\link{.deprecate}()}
#'
#' @noRd
#'
.deprecation_stage <- function(when) {
  stopifnot(rlang::is_character(when, n = 1L) && nzchar(x = when))
  when <- numeric_version(when, strict = FALSE)
  if (rlang::is_na(when)) {
    rlang::abort("'when' must be a valid version")
  }
  # Get list of releases
  releases <- as.data.frame(read.dcf(system.file(
    "extdata",
    "releases.dcf",
    package = .pkgenv$pkgname,
    lib.loc = .pkgenv$libname,
    mustWork = TRUE
  )))
  # Add fake release for current minor release if needed
  mm <- unique(vapply(
    X = releases$Version,
    FUN = function(x) {
      x <- unlist(strsplit(x = x, split = "\\."))
      return(paste0(x[1:2], collapse = "."))
    },
    FUN.VALUE = character(length = 1L),
    USE.NAMES = FALSE
  ))
  cmm <- unlist(strsplit(
    x = as.character(utils::packageVersion(.pkgenv$pkgname)),
    split = "\\."
  ))
  if (all(paste0(cmm[1:2], collapse = ".") > mm)) {
    dates <- file.info(list.files(
      base::system.file(package = .pkgenv$pkgname),
      full.names = TRUE,
      recursive = TRUE
    ))$mtime
    releases <- rbind(
      releases,
      data.frame(
        Version = sprintf(fmt = "%s.%s.0", cmm[1L], cmm[2L]),
        Date = format(as.POSIXlt(dates[which.max(dates)]), format = "%Y-%m-%d")
      )
    )
  }
  # Check to see if the deprecation is scheduled for a future release
  # If so, exit out
  if (all(when > releases$Version)) {
    return(invisible(NULL))
  }
  if (!when %in% releases$Version) {
    stop(sprintf(
      fmt = "Unknown %s release: '%s'",
      .pkgenv$pkgname,
      as.character(when)
    ))
  }
  # Get current package version
  current <- .tiledbsoma_deprecation_version()
  if (current < when) {
    # Deprecation will happen in the future
    return(invisible(NULL))
  }
  # Get current package's release date
  if (current %in% releases$Version) {
    date <- as.POSIXlt(releases[releases$Version == current, "Date"])
  } else {
    date <- as.POSIXlt(read.dcf(
      system.file(
        "DESCRIPTION",
        package = .pkgenv$pkgname,
        lib.loc = .pkgenv$libname,
        mustWork = TRUE
      ),
      fields = "Date/Publication"
    ))
    if (is.na(date)) {
      dates <- file.info(list.files(
        base::system.file(package = .pkgenv$pkgname),
        full.names = TRUE,
        recursive = TRUE
      ))$mtime
      date <- as.POSIXlt(dates[which.max(dates)])
    }
  }
  weeks <- as.double(
    difftime(
      time1 = date,
      time2 = as.POSIXlt(x = releases[releases$Version == when, "Date"]),
      units = "weeks"
    ),
    units = "weeks"
  )
  .as_integer_version <- function(x) {
    x <- vapply(
      X = unlist(strsplit(as.character(x), split = "\\.")),
      FUN = as.integer,
      FUN.VALUE = integer(length = 1L)
    )
    if (length(x) > 4L) {
      x <- x[1:4]
    }
    names(x) <- c("major", "minor", "patch", "devel")[seq_along(x)]
    return(x)
  }
  cc <- .as_integer_version(current)
  ww <- .as_integer_version(when)
  defunct <- cc[["major"]] > ww[["major"]] || (
    weeks > 12L &&
      cc[["major"]] == ww[["major"]] && (cc[["minor"]] - ww[["minor"]]) >= 2L
  )
  if (defunct) {
    return("defunct")
  }
  return("deprecate")
}
