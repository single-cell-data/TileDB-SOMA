
string_collapse <- function(x, sep = ", ") {
  return(glue::glue_collapse(x, sep = sep, width = getOption("width", Inf)))
}

n_unique <- function(x) {
  length(unique(x))
}

vapply_char <- function(X, FUN, ..., USE.NAMES = TRUE) {
  vapply(X, FUN, FUN.VALUE = character(1L), ..., USE.NAMES = USE.NAMES)
}

vapply_lgl <- function(X, FUN, ..., USE.NAMES = TRUE) {
  vapply(X, FUN, FUN.VALUE = logical(1L), ..., USE.NAMES = USE.NAMES)
}

vapply_int <- function(X, FUN, ..., USE.NAMES = TRUE) {
  vapply(X, FUN, FUN.VALUE = integer(1L), ..., USE.NAMES = USE.NAMES)
}

rename <- function(x, names) {
  stopifnot(
    "'x' must be named" = is_named(x),
    "'names' must be a named character vector" = is_named(names),
    "All 'names' must be in 'x'" = all(names %in% names(x))
  )

  name_index <- match(names, names(x))
  names(x)[name_index] <- names(names)
  x
}

# Return y if x is NULL, else x
`%||%` <- function(x, y) {
  if (missing(x) || is.null(x) || length(x) == 0) y else x
}

err_to_warn <- function(err, immediate. = TRUE) {
  warning(conditionMessage(err), call. = FALSE, immediate. = immediate.)
  return(invisible(err))
}

null <- function(...) {
  return(NULL)
}

random_name <- function(length = 5L, chars = letters, ...) {
  stopifnot(
    "'length' must be a single integer" = rlang::is_integerish(length, n = 1L),
    "'chars' must be character" = is.character(chars)
  )
  chars <- unique(unlist(strsplit(chars, split = "")))
  return(paste(sample(chars, size = length, ...), collapse = ""))
}

uns_hint <- function(type = c("1d", "2d")) {
  type <- match.arg(type)
  hint <- list(paste0("array_", type))
  names(hint) <- "soma_uns_outgest_hint"
  return(hint)
}

.encode_as_char <- function(x) {
  return(switch(
    EXPR = typeof(x),
    double = sprintf("%a", x),
    x
  ))
}

.err_to_warn <- function(err) {
  warning(warningCondition(
    message = conditionMessage(err),
    class = setdiff(class(err), c("warning", "simpleError", "error", "condition")),
    call = conditionCall(err)
  ))
}

.decode_from_char <- function(x) {
  stopifnot(is.character(x))
  double <- paste0(
    "^",
    c(
      "[-]?0x[0-9a-f](\\.[0-9a-f]+)?p[+-][0-9]+",
      "[-]?Inf",
      "NA",
      "NaN"
    ),
    "$",
    collapse = "|"
  )
  return(if (all(grepl(double, x))) {
    as.numeric(x)
  } else {
    x
  })
}

#' Is an Object Integerish
#'
#' @inheritParams rlang::is_integerish
#'
#' @return \code{TRUE} if \code{x} is integerish, otherwise \code{FALSE}
#'
#' @keywords internal
#'
#' @noRd
#'
.is_integerish <- function(x, n = NULL, finite = NULL) {
  UseMethod(generic = ".is_integerish", object = x)
}

#' @method .is_integerish default
#' @export
#'
.is_integerish.default <- function(x, n = NULL, finite = NULL) {
  return(rlang::is_integerish(x = x, n = n, finite = finite))
}

#' @method .is_integerish integer64
#' @export
#'
.is_integerish.integer64 <- function(x, n = NULL, finite = NULL) {
  res <- if (!is.null(x = n)) {
    stopifnot(
      "'n' must be a single integerish value" = .is_integerish(x = n) &&
        length(x = n) == 1L &&
        is.finite(x = n)
    )
    length(x = x) == n
  } else {
    TRUE
  }
  res <- res && if (!is.null(x = finite)) {
    stopifnot(isTRUE(x = finite) || isFALSE(x = finite))
    # In `rlang::is_integerish()`,
    # `finite = TRUE`: all values are finite
    # `finite = FALSE`: at least one value is infinite
    # `bit64::is.infinite()` returns FALSE for NA
    ifelse(
      test = finite,
      yes = all(is.finite(x = x)),
      no = any(is.infinite(x = x) | is.na(x = x))
    )
  } else {
    TRUE
  }
  return(res)
}

#' @method .is_integerish Field
#' @export
#'
.is_integerish.Field <- function(x, n = NULL, finite = NULL) {
  return(.is_integerish(x = x$type, n = n, finite = finite))
}

#' @method .is_integerish Array
#' @export
#'
.is_integerish.Array <- .is_integerish.Field

#' @method .is_integerish ChunkedArray
#' @export
#'
.is_integerish.ChunkedArray <- .is_integerish.Field

#' @method .is_integerish DataType
#' @export
#'
.is_integerish.DataType <- function(x, n = NULL, finite = NULL) {
  return(grepl(pattern = "^[u]?int[[:digit:]]{1,2}$", x = x$name))
}

.maybe_muffle <- function(w, cond = getOption("verbose", default = FALSE)) {
  if (isTRUE(x = cond)) {
    warning(warningCondition(
      message = conditionMessage(w),
      class = setdiff(class(w), c("warning", "simpleError", "error", "condition")),
      call = conditionCall(w)
    ))
  } else {
    tryInvokeRestart("muffleWarning")
  }
}

#' Generate a SOMA Metadata Type Hint
#'
#' @param type A character vector giving the class of an object
#'
#' @return A named list where the name is \dQuote{\code{soma_r_type_hint}} and
#' the value is a single string value giving the R class; this value changes
#' based on the type of object a type hint is being generated for
#' \describe{
#'  \item{Simple S3 objects}{the R class (eg. \dQuote{\code{data.frame}})}
#'  \item{S3 classes with inheritance}{a JSON-array encoding of the R class
#'   (eg. \dQuote{\code{["matrix", "array"]}})}
#'  \item{S4 classes}{the R package and class encoded in
#'   \dQuote{\code{pkg:class}} form (eg. \dQuote{\code{Matrix:dgCMatrix}})}
#' }
#'
#' @keywords internal
#'
#' @examples
#' # Type hint for S3 classes
#' .type_hint("data.frame") # data.frame
#'
#' # Type hint for complex S3 classes
#' .type_hint(class(matrix())) # ["matrix","array"]
#'
#' # Type hint for S4 classes
#' .type_hint(class(Matrix::Matrix())) # Matrix::ldiMatrix
#'
#' @noRd
#'
.type_hint <- function(type) {
  nm <- 'soma_r_type_hint'
  if (is.null(type)) {
    hint <- list(NULL)
    names(hint) <- nm
    return(hint)
  }
  stopifnot(
    "'type' must be a non-empty character" = is.character(type) && all(nzchar(type))
  )
  def <- if (length(type) > 1L) {
    paste0('[', paste(dQuote(type, FALSE), collapse = ','), ']')
  } else {
    tryCatch(
      expr = methods::getClassDef(type),
      error = function(e) type
    )
  }
  if (inherits(def, c('classUnionRepresentation', 'refClassRepresentation'))) {
    def <- sprintf('%s:%s', def@package, def@className)
  } else if (inherits(def, 'classRepresentation')) {
    btypes <- vapply_char(
      X = gsub(
        pattern = '^is\\.',
        replacement = '',
        x = grep(
          pattern = '<-',
          x = grep(
            pattern = '^is\\.',
            x = utils::lsf.str(envir = baseenv()),
            value = TRUE
          ),
          value = TRUE,
          invert = TRUE
        )
      ),
      FUN = function(x) ifelse(
        test = grepl(pattern = '^data\\.frame', x = x),
        yes = paste(strsplit(x, split = '\\.')[[1L]][1:2], collapse = '.'),
        no = strsplit(x, split = '\\.')[[1L]][1L]
      ),
      USE.NAMES = FALSE
    )
    def <- switch(
      EXPR = def@package,
      methods = {
        def <- if ('oldClass' %in% names(def@contains) || def@className %in% btypes) {
          as.character(def@className)
        } else {
          sprintf('%s:%s', def@package, def@className)
        }
        def
      },
      sprintf(fmt = '%s:%s', def@package, def@className)
    )
  }
  return(list(soma_r_type_hint = def))
}

#' Read the SOMA Join IDs from an Array
#'
#' @param x A \code{\link{SOMASparseNDarray}} or \code{\link{SOMADataFrame}}
#'
#' @return An \code{\link[bit64]{integer64}} vector with the SOMA join IDs
#'
#' @keywords internal
#'
#' @noRd
#'
.read_soma_joinids <- function(x, ...) {
  stopifnot(inherits(x = x, what = "SOMAArrayBase"))
  oldmode <- x$mode()
  on.exit(
    x$reopen(oldmode, tiledb_timestamp = x$tiledb_timestamp),
    add = TRUE,
    after = FALSE
  )
  op <- options(arrow.int64_downcast = FALSE)
  on.exit(options(op), add = TRUE, after = FALSE)
  ids <- UseMethod(generic = ".read_soma_joinids", object = x)
  return(ids)
}

#' @rdname dot-read_soma_joinids
#'
#' @noRd
#'
#' @method .read_soma_joinids SOMADataFrame
#' @export
#'
.read_soma_joinids.SOMADataFrame <- function(x, ...) {
  x$reopen("READ", tiledb_timestamp = x$tiledb_timestamp)
  return(x$read(column_names = "soma_joinid")$concat()$GetColumnByName("soma_joinid")$as_vector())
}

#' @param axis Which dimension to read (zero-based)
#'
#' @rdname dot-read_soma_joinids
#'
#' @noRd
#'
#' @method .read_soma_joinids SOMASparseNDArray
#' @export
#'
.read_soma_joinids.SOMASparseNDArray <- function(x, axis = 0L, ...) {
  stopifnot(
    "'axis' must be a single positive integer" = is.integer(axis) &&
      length(axis) == 1L
  )
  if (axis < 0L || axis >= length(x$dimnames())) {
    stop("'axis' must be between 0 and ", length(x$dimnames()), call. = FALSE)
  }
  x$reopen("READ", tiledb_timestamp = x$tiledb_timestamp)
  dimname <- x$dimnames()[axis + 1L]
  sr <- mq_setup(
    uri = x$uri,
    soma_context(),
    colnames = dimname,
    timestamprange = x$.tiledb_timestamp_range
  )
  return(TableReadIter$new(sr)$concat()$GetColumnByName(dimname)$as_vector())
}

#' Pad Names of a Character Vector
#'
#' Fill in missing names of a vector using missing values of said vector
#'
#' @param x A character vector
#'
#' @return \code{x} with any missing names set to the values of \code{x}; see
#' examples for more details
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' x1 <- c("a", "b", "c")
#' pad_names(x1) # returns c(a = "a", b = "b", c = "c")
#'
#' x2 <- c(a = "x", b = "y", c = "z")
#' pad_names(x2) # returns c(a = "x", b = "y", c = "z")
#'
#' x3 <- c(a = "x", "y", c = "z")
#' pad_names(x3) # returns c(a = "x", y = "y", c = "z")
#'
pad_names <- function(x) {
  stopifnot(
    is.character(x)
  )
  if (is.null(names(x))) {
    return(stats::setNames(nm = x))
  }
  unnamed <- !nzchar(names(x))
  names(x)[unnamed] <- x[unnamed]
  return(x)
}

# For use in read-only R6 active bindings
read_only_error <- function(field_name) {
  stop(
    sprintf("'%s' is a read-only field.", field_name),
    call. = FALSE
  )
}

SOMA_OBJECT_TYPE_METADATA_KEY <- "soma_object_type"
SOMA_ENCODING_VERSION_METADATA_KEY <- "soma_encoding_version"
SOMA_ENCODING_VERSION <- "1.1.0"

# This is for internal logging purposes. Context:
# * We have R (and Python) code with function names the user invokes.
# * These call C++ functions which can throw their own error messages.
# * It's crucial that the C++ code "knows" the name of the function
#   as typed by the user, not whatever (possibly out-of-date) guess
#   the C++ code may have.
.name_of_function <- function() {
  # Tricky bits:
  # * This might be `obj$foo`
  # * The sys.call can return a parse-tree component (typeof = language)
  #   with the '$', 'obj', and 'foo' -- hence the as.character
  # * Even then there can be a second component returned like 'c(1)'
  #   -- hence the [[1]]
  # * Then remove the 'obj$' from 'obj$foo'
  name <- as.character(sys.call(sys.parent(n=1)))[[1]]
  name <- sub('.*\\$', replacement = '', x = name)
  return(name)
}
