
##' Convenience Wrappers for 'RcppSpdlog' Logging From 'spdlog'
##'
##' Several short wrappers for functions from 'RcppSpdlog' package are provided
##' as a convenience.  Given the potential for clashing names of common and popular
##' functions names we do \emph{not} recommend the import the whole package but rather
##' do \code{importFrom(RcppSpdlog, set_pattern)} (or maybe \code{importFrom(RcppSpdlog,
##' set_pattern)}). After that, functionality can be accessed via a convenient shorter
##' form such as for example \code{spld::info()} to log at the \sQuote{info} level.
##'
##' @param name Character value for the name of the logger instance
##' @param level Character value for the logging level
##' @param s Character value for pattern, level, or logging message
##' @return Nothing is returned from these functions as they are invoked for their side-effects.
##' @examples
##' spdl_setup("exampleDemo", "warn")
##' spdl_info("Not seen as level 'info' below 'warn'")
##' spdl_warn("This warning message is seen")
##' spdl_set_level("info")
##' spdl_info("Now this informational message is seen too")
##' spdl_info(sprintf("Calls to sprintf() can be used %s as we see %d", "as well", 42L))
spdl_setup       <- function(name = "default", level = "warn") RcppSpdlog::log_setup(name, level)

##' @rdname spdl_setup
spdl_drop        <- function(s)  RcppSpdlog::log_drop(s)

##' @rdname spdl_setup
spdl_set_pattern <- function(s)  RcppSpdlog::log_set_pattern(s)

##' @rdname spdl_setup
spdl_set_level   <- function(s)  RcppSpdlog::log_set_level(s)

##' @rdname spdl_setup
spdl_trace       <- function(s)  RcppSpdlog::log_trace(s)

##' @rdname spdl_setup
spdl_debug       <- function(s)  RcppSpdlog::log_debug(s)

##' @rdname spdl_setup
spdl_info        <- function(s)  RcppSpdlog::log_info(s)

##' @rdname spdl_setup
spdl_warn        <- function(s)  RcppSpdlog::log_warn(s)

##' @rdname spdl_setup
spdl_error       <- function(s)  RcppSpdlog::log_error(s)

##' @rdname spdl_setup
spdl_critical    <- function(s)  RcppSpdlog::log_critical(s)
