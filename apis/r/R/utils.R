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

get_soma_context <- function(context, tiledbsoma_ctx, what = NULL) {
  if (!is.null(tiledbsoma_ctx)) {
    .deprecate(
        what=what,
        when="2.3.0",
        details="Use `context` instead."
    )
  }
  if (is.null(context)) {
    if (is.null(tiledbsoma_ctx)) {
        context <- .pkgenv[["somactx"]]
        if (is.null(context)) {
          return(set_default_context())
        }
        reurn(context)
    }
    if (!inherits(x = tiledbsoma_ctx, what = 'SOMATileDBContext')) {
      stop(
        "'tiledbsoma_ctx' must be a SOMATileDBContext object",
        call. = FALSE
      )
    }
    return(SOMAContext$new(config = unlist(tiledbsoma_ctx$to_list())))
  }
  if (!is.null(tiledbsoma_ctx)) {
    warning(
      "Both 'context' and 'tiledbsoma_ctx' were provided, using 'context' only"
    )
  }
  if (!inherits(x = context, what = 'SOMAContext')) {
    stop(
      "'context' must be a SOMAContext object",
      call. = FALSE
    )
  }
  return(context)
}


#' Generate a Block Size for Matrix Iteration
#'
#' Generate block sizes for matrix iteration; the block sizes are calculated
#' while keeping one axis static and chunking on the alternate axis
#'
#' \code{context} is an option available to \code{write_soma()}; we can
#' use the options in the context to determine how much memory we can use from
#' the option \dQuote{\code{soma.init_buffer_bytes}}. A default value of
#' \eqn{33,554,432} bytes is used when no context is provided, or the context
#' does not specify a memory size. This value is taken from the CI workflows
#' for \pkg{tiledbsoma}
#'
#' @param n Number of entries on the static (non-iterated) axis
#' @param context A \code{SOMAContext} object
#'
#' @return Number of entries on the alternate axis to chunk
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examplesIf requireNamespace("SeuratObject", quietly = TRUE)
#' data("pbmc_small", package = "SeuratObject")
#' pbmc_small <- UpdateSeuratObject(pbmc_small)
#'
#' # Pull a matrix from `pbmc_small`, transpose to be obs x var
#' # instead of var x obs
#' mat <- t(pbmc_small[["RNA"]]$counts)
#' n_var <- ncol(mat)
#'
#' # Use a context to limit tp half a megabyte
#' ctx <- SOMAContext$new(c(
#'   soma.init_buffer_bytes = as.character(0.5 * (1024L ^ 2L))
#' ))
#'
#' # Generate a block size to iterate across the obs axis
#' # `n` is the number of entries on the static axis, not the iterated axis
#' .block_size(n = n_var, context = ctx)
#'
.block_size <- function(n, context) {
  if (!rlang::is_integerish(n, n = 1L, finite = TRUE) || n <= 0L) {
    rlang::abort("'n' must be a single, finite, positive integer value")
  }
  if (!inherits(context, "SOMAContext")) {
    rlang::abort("'context' must be a SOMAContext object")
  }

  # Try to get the "soma.init_buffer_bytes option from the context.
  # - If it wasn't set, then use a default value.
  # - If it was set but isn't numeric, then throw an error.
  bytes <- 33554432 # default value - over-write if value in config
  config <- context$get_config()
  bytes_key = "soma.init_buffer_bytes"
  if (bytes_key %in% names(config)) {
    bytes <- tryCatch(expr = as.numeric(config[bytes_key], warning = stop))
  }

  # Calculate the blocks size assuming we're working with numeric matrices
  # Not integer or logical
  num_bytes <- utils::object.size(numeric(1L))
  return(floor(bytes / as.numeric(num_bytes) / n))
}

.encode_as_char <- function(x) {
  return(switch(EXPR = typeof(x), double = sprintf("%a", x), x))
}

.err_to_warn <- function(err) {
  warning(warningCondition(
    message = conditionMessage(err),
    class = setdiff(
      class(err),
      c("warning", "simpleError", "error", "condition")
    ),
    call = conditionCall(err)
  ))
}

.decode_from_char <- function(x) {
  stopifnot(is.character(x))
  double <- paste0(
    "^",
    c("[-]?0x[0-9a-f](\\.[0-9a-f]+)?p[+-][0-9]+", "[-]?Inf", "NA", "NaN"),
    "$",
    collapse = "|"
  )
  return(
    if (all(grepl(double, x))) {
      as.numeric(x)
    } else {
      x
    }
  )
}

#' Is an Object a Domain Specification
#'
#' Check that an object is a domain specification. Valid domain specifications
#' are named lists where each entry is either \code{NULL} or a two-length
#' atomic vector (no lists or factors) giving the bounds
#'
#' @param x An \R object
#' @param nm Vector of expected names of the domain specification
#'
#' @return \code{TRUE} if \code{x} is a valid domain specification,
#' otherwise \code{FALSE}
#'
#' @keywords internal
#'
#' @noRd
#'
.is_domain <- function(x, nm) {
  if (!is.character(nm)) {
    stop("'nm' must be a character vector")
  }
  return(
    rlang::is_list(x = x, n = length(x = nm)) &&
      identical(x = sort(names(x = x)), y = sort(nm)) &&
      all(vapply(
        X = x,
        FUN = function(i) {
          return(
            is.null(i) ||
              (is.atomic(i) && !is.factor(i) && length(i) == 2L)
          )
        },
        FUN.VALUE = logical(1L)
      ))
  )
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
  res <- res &&
    if (!is.null(x = finite)) {
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
      class = setdiff(
        class(w),
        c("warning", "simpleError", "error", "condition")
      ),
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
    "'type' must be a non-empty character" = is.character(type) &&
      all(nzchar(type))
  )
  def <- if (length(type) > 1L) {
    paste0('[', paste(dQuote(type, FALSE), collapse = ','), ']')
  } else {
    tryCatch(expr = methods::getClassDef(type), error = function(e) type)
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
      FUN = function(x) {
        ifelse(
          test = grepl(pattern = '^data\\.frame', x = x),
          yes = paste(strsplit(x, split = '\\.')[[1L]][1:2], collapse = '.'),
          no = strsplit(x, split = '\\.')[[1L]][1L]
        )
      },
      USE.NAMES = FALSE
    )
    def <- switch(
      EXPR = def@package,
      methods = {
        def <- if (
          'oldClass' %in% names(def@contains) || def@className %in% btypes
        ) {
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
  return(x$read(column_names = "soma_joinid")$concat()$GetColumnByName(
    "soma_joinid"
  )$as_vector())
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
    x$context$handle,
    colnames = dimname,
    timestamprange = x$.tiledb_timestamp_range
  )
  return(TableReadIter$new(sr)$concat()$GetColumnByName(dimname)$as_vector())
}

#' Determine if an S3 Method Exists
#'
#' Check to see that S3 method dispatch could reasonably be performed. Exposed
#' S3 method selectors do not follow S4 inheritance, so asking if a method
#' exists using \code{\link[utils:getS3method]{utils::getS3method}()} will not
#' always work. This function simply checks to see if and S3 method
#' \emph{could} be selected, but will not attempt to select one for use
#'
#' @param f Name of generic function
#' @param class Name of class
#'
#' @return \code{TRUE} if an S3 method \code{<f>.<class>()} can be found,
#' otherwise \code{FALSE}
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' utils::methods("write_soma") # list all methods for `write_soma()`
#'
#' # `utils::.getS3method()` reports that there is no method for dgCMatrices
#' # even though `write_soma.Matrix()` is defined
#' tryCatch(utils::getS3method("write_soma", "dgCMatrix"), error = force)
#'
#' # `.s3_method_defined()` will report that a method for dgCMatrices exists
#' .s3_method_defined("write_soma", "dgCMatrix") # TRUE
#'
#' # For classes where there is no method, `.s3_method_defined()` returns FALSE
#' .s3_method_defined("write_soma", "classRepresentation")
#' .s3_method_defined("write_soma", "vector")
#'
.s3_method_defined <- function(f, class) {
  stopifnot(
    rlang::is_character(f, n = 1L) && nzchar(f),
    rlang::is_character(class) && all(nzchar(class))
  )
  methods <- attr(utils::methods(generic.function = f), which = "info")
  if (!nrow(methods)) {
    rlang::abort(sprintf("'%s' is not a generic function", f))
  }
  methods$method <- row.names(methods)
  methods <- methods[!methods$isS4, , drop = FALSE]
  if (!nrow(methods)) {
    rlang::abort(sprintf("'%s' is not an S3 generic", f))
  }
  methods <- do.call(
    what = rbind,
    args = lapply(
      X = split(x = methods, f = methods$from),
      FUN = function(df) {
        fakes <- tools::nonS3methods(unique(df$from))
        return(df[!df$method %in% fakes, , drop = FALSE])
      }
    )
  )
  if (!nrow(methods)) {
    rlang::abort(sprintf("No actual methods found for '%s'", f))
  }
  row.names(methods) <- NULL
  if (length(class) == 1L) {
    cdef <- methods::getClassDef(Class = class)
    if (!is.null(cdef) && !"oldClass" %in% names(cdef@contains)) {
      class <- c(class, names(cdef@contains))
    }
  }
  return(any(sprintf("%s.%s", f, class) %in% methods$method))
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
  stopifnot(is.character(x))
  if (is.null(names(x))) {
    return(stats::setNames(nm = x))
  }
  unnamed <- !nzchar(names(x))
  names(x)[unnamed] <- x[unnamed]
  return(x)
}

# For use in read-only R6 active bindings
read_only_error <- function(field_name) {
  stop(sprintf("'%s' is a read-only field.", field_name), call. = FALSE)
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
  name <- as.character(sys.call(sys.parent(n = 1)))[[1]]
  name <- sub('.*\\$', replacement = '', x = name)
  return(name)
}
