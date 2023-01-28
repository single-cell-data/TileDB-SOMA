#' @include MappingBase.R
#'
NULL

PlatformConfig <- R6::R6Class(
  classname = 'PlatformConfig',
  inherit = MappingBase,
  public = list(
    initialize = function() {
      .NotYetImplemented()
    },
    platforms = function() {
      return(self$keys())
    },
    ops = function(platform = NULL) {
      platform <- platform %||% self$platforms()[1L]
      if (isTRUE(x = platform)) {
        ops <- Reduce(
          f = union,
          x = lapply(
            X = self$platforms(),
            FUN = \(p) private$.data[[p]]$keys()
          )
        )
        return(ops)
      }
      platform <- platform[1L]
      platform <- match.arg(arg = platform, choices = self$platforms())
      return(super$get(key = platform)$keys())
    },
    get_ops = function(platform) {
      platform <- platform[1L] %||% self$platforms()[1L]
      platform <- match.arg(arg = platform, choices = self$platforms())
      return(super$get(key = platform))
    },
    set = function(platform, op, value) {
      stopifnot(
        is.character(x = platform) && length(x = platform) == 1L,
        is.character(x = op) && length(x = op) == 1L,
        inherits(x = value, what = 'ScalarMap')
      )
      platform <- tryCatch(
        expr = match.arg(arg = platform, choices = self$platform()),
        error = \(...) platform
      )
      pconf <- if (platform %in% self$.platforms()) {
        self$get(key = platform)
      } else {
        ConfigList$new()
      }
      pconf$add_map(key = op, value = value)
      return(invisible(x = self))
    }
  ),
  private = list()
)

# setValidity(
#   Class = 'PlatformConfig',
#   method = function(object) {
#     valid <- NULL
#     if (!is_named_list(x = object)) {
#       valid <- c(valid, "All entries must be named")
#     }
#     # Iterate through platforms
#     for (i in seq_along(along.with = object)) {
#       ii <- names(x = object)[i]
#       ni <- if (isTRUE(x = nzchar(x = ii))) {
#         sQuote(x = ii)
#       } else {
#         paste('#', i)
#       }
#       val <- object[[i]]
#       if (!is_named_list(x = val)) {
#         valid <- c(
#           valid,
#           paste0(
#             "All top-level entries must be named lists (offending entry: ",
#             ni,
#             ")"
#           )
#         )
#       }
#       if (!is.list(x = val)) {
#         next
#       }
#       # Iterate through operations
#       for (j in seq_along(along.with = val)) {
#         jj <- names(x = val)[j]
#         nj <- if (isTRUE(x = nzchar(x = jj))) {
#           sQuote(x = jj)
#         } else {
#           paste('#', j)
#         }
#         vj <- val[[j]]
#         if (!length(x = vj)) {
#           next
#         }
#         if (!is_named_list(x = vj)) {
#           valid <- c(
#             valid,
#             paste0(
#               "All entries in ",
#               ni,
#               " must be named lists (offending entry: ",
#               nj,
#               ")"
#             )
#           )
#         }
#         if (!is.list(x = vj)) {
#           next
#         }
#         # Iterate through options
#         opt_check <- vapply(
#           X = vj,
#           FUN = \(x) length(x = x) == 1L && is.atomic(x = x),
#           FUN.VALUE = logical(length = 1L)
#         )
#         if (!all(opt_check)) {
#           valid <- c(
#             valid,
#             paste0(
#               "All entries in ",
#               ni,
#               ':',
#               nj,
#               " must be single atomic values"
#             )
#           )
#         }
#       }
#     }
#     return(valid %||% TRUE)
#   }
# )
