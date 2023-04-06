#' Write a SOMA Object from an \R Object
#'
#' Convert \R objects to their appropriate SOMA counterpart
#' function and methods can be written for it to provide a high-level
#' \R \eqn{\rightarrow} SOMA interface
#'
#' @param x An object
#' @param uri URI for resulting SOMA object
#' @template param-dots-method
#' @param platform_config Optional \link[tiledbsoma:PlatformConfig]{platform
#' configuration}
#' @param tiledbsoma_ctx Optional \code{\link{SOMATileDBContext}}
#'
#' @return A \code{\link{SOMAExperiment}} with the data from \code{x}
#'
#' @section Known methods:
#' \itemize{
#'  \item \link[tiledbsoma:write_soma.Seurat]{Writing Seurat objects}
#' }
#'
#' @export
#'
#' @noMd
#'
write_soma <- function(x, uri, ..., platform_config = NULL, tiledbsoma_ctx = NULL) {
  UseMethod(generic = 'write_soma', object = x)
}

#' Write R Objects to SOMA
#'
#' Various helpers to write R objects to SOMA
#'
#' @inheritParams write_soma
#' @param soma_parent The parent \link[tiledbsoma:SOMACollection]{collection} (eg. a
#' \code{\link{SOMACollection}}, \code{\link{SOMAExperiment}}, or
#' \code{\link{SOMAMeasurement}})
#' @param absolute \strong{\[Internal use only\]} Is \code{uri} absolute
#'
#' @return The resulting SOMA \link[tiledbsoma:SOMASparseNDArray]{array} or
#' \link[tiledbsoma:SOMADataFrame]{data frame}
#'
#' @name write_soma_objects
#' @rdname write_soma_objects
#'
#' @keywords internal
#'
NULL

#' @param index_name The name of the column in \code{x} with the index (row names);
#' by default, will automatically add the row names of \code{x} to an attribute
#' named \dQuote{\code{index}} to the resulting \code{\link{SOMADataFrame}}
#'
#' @rdname write_soma_objects
#'
#' @section Writing Data Frames:
#' \link[base:data.frame]{Data frames} are written out as
#' \code{\link{SOMADataFrame}s}. The following transformations
#' are applied to \code{x}:
#' \itemize{
#'  \item row names are added to a column in \code{x} entitled
#'   \dQuote{\code{index}}, \dQuote{\code{_index}}, or a random name if
#'   either option is already present in \code{x}
#'  \item a column \dQuote{\code{soma_joinid}} will be automatically added
#'   going from \code{[0, nrow(x) - 1]} encoded as
#'   \link[bit64:integer64]{64-bit integers}; any columns called
#'   \dQuote{\code{soma_joinid}} will be overwritten if present in \code{x}
#'  \item all factor columns will be coerced to characters
#'  \item all columns not containing \link[base:is.atomic]{atomic} types
#'   (excluding \link[base:factor]{factors}, \code{\link[base]{complex}es},
#'   and \code{\link[base]{raw}s}) will be removed
#' }
#' The array type for each column will be determined by
#' \code{\link[arrow:infer_type]{arrow::infer_type}()}
#'
#' @method write_soma data.frame
#' @export
#'
write_soma.data.frame <- function(
  x,
  uri,
  soma_parent,
  index_name = NULL,
  ...,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  absolute = FALSE
) {
  stopifnot(
    "'index_name' must be a single character value" = is.null(index_name) ||
      is_scalar_character(index_name)
  )
  # Create a proper URI
  uri <- .check_soma_uri(
    uri = uri,
    soma_parent = soma_parent,
    absolute = absolute
  )
  # Clean up data types in `x`
  remove <- vector(mode = 'logical', length = ncol(x))
  for (i in seq_len(ncol(x))) {
    col <- names(x)[i]
    if (is.factor(x[[col]])) {
      x[[col]] <- as.character(x[[col]])
    }
    remove[i] <- !inherits(
      x = try(expr = arrow::infer_type(x[[col]]), silent = TRUE),
      what = 'DataType'
    )
    # remove[i] <- !inherits(
    #   x = x[[col]],
    #   what = setdiff(x = .SCALAR_TYPES(), y = 'any')
    # )
  }
  if (all(remove)) {
    stop(
      "None of the columns in the data frame can be written out",
      call. = FALSE
    )
  } else if (any(remove)) {
    warning(
      paste(
        strwrap(paste(
          "Removing the following columns due to incompatible data type:",
          paste(sQuote(names(x)[remove]), collapse = ', ')
        )),
        collapse = '\n'
      ),
      call. = FALSE,
      immediate. = TRUE
    )
    x <- x[, !remove, drop = FALSE]
  }
  # Check `index_name`
  index_name <- index_name %||% attr(x = x, which = 'index')
  if (is.null(index_name)) {
    x <- .df_index(x = x, ...)
    index_name <- attr(x = x, which = 'index')
  }
  if (!index_name %in% names(x)) {
    stop(
      "Unable to find ",
      sQuote(index_name),
      " in the provided data frame",
      call. = FALSE
    )
  }
  # Add `soma_joinid` to `x`
  if ('soma_joinid' %in% names(x)) {
    warning(
      "Found column 'soma_joinid' in the data frame, overwriting with new 'soma_joinid'",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  x$soma_joinid <- bit64::seq.integer64(from = 0L, to = nrow(x) - 1L)
  # Create the SOMADataFrame
  tbl <- arrow::arrow_table(x)
  sdf <- SOMADataFrameCreate(
    uri = uri,
    schema = tbl$schema,
    index_column_names = 'soma_joinid',
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  # Write and return
  sdf$write(tbl)
  return(sdf)
}

#' @param sparse Create a \link[tiledbsoma:SOMASparseNDArray]{sparse} or
#' \link[tiledbsoma:SOMADenseNDArray]{dense} array from \code{x}
#' @param type \link[arrow:data-type]{Arrow type} for encoding \code{x}
#' (eg. \code{\link[arrow:data-type]{arrow::int32}()}); by default, attempts to
#' determine arrow type with \code{\link[arrow:infer_type]{arrow::infer_type}()}
#' @param transpose Transpose \code{x} before writing
#'
#' @rdname write_soma_objects
#'
#' @section Writing Dense Matrices:
#' Dense matrices are written as two-dimensional
#' \link[tiledbsoma:SOMADenseNDArray]{dense arrays}. The overall shape of the
#' array is determined by \code{dim(x)} and the type of the array is determined
#' by \code{type} or \code{\link[arrow:infer_type]{arrow::infer_type}(x)}
#'
#' @method write_soma matrix
#' @export
#'
write_soma.matrix <- function(
  x,
  uri,
  soma_parent,
  sparse = TRUE,
  type = NULL,
  transpose = FALSE,
  ...,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  absolute = FALSE
) {
  stopifnot(
    "'sparse' must be a single logical value" = is_scalar_logical(sparse),
    "'type' must be an Arrow type" = is.null(type) ||
      (R6::is.R6(type) && inherits(x = type, what = 'DataType')),
    "'transpose' must be a single logical value" = is_scalar_logical(transpose)
  )
  if (!isTRUE(sparse) && inherits(x = x, what = 'sparseMatrix')) {
    msg <- "A sparse matrix was provided and a dense array was asked for, creating a sparse array"
    warning(paste(strwrap(msg), collapse = '\n'), call. = FALSE, immediate. = TRUE)
    sparse <- TRUE
  }
  # Create a sparse array
  if (isTRUE(sparse)) {
    return(write_soma(
      x = methods::as(object = x, Class = 'TsparseMatrix'),
      uri = uri,
      soma_parent = soma_parent,
      type = type,
      transpose = transpose,
      ...,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      absolute = absolute
    ))
  }
  # Create a dense array
  if (inherits(x = x, what = 'Matrix')) {
    x <- as.matrix(x)
  }
  # Create a proper URI
  uri <- .check_soma_uri(
    uri = uri,
    soma_parent = soma_parent,
    absolute = absolute
  )
  # Transpose the matrix
  if (isTRUE(transpose)) {
    x <- t(x)
  }
  # Create the array
  array <- SOMADenseNDArrayCreate(
    uri = uri,
    type = type %||% arrow::infer_type(x),
    shape = dim(x),
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  # Write and return
  array$write(x)
  return(array)
}

#' @rdname write_soma_objects
#'
#' @method write_soma Matrix
#' @export
#'
write_soma.Matrix <- write_soma.matrix

#' @rdname write_soma_objects
#'
#' @section Writing Sparse Matrices:
#' Sparse matrices are written out as two-dimensional
#' \link[tiledbsoma:SOMASparseNDArray]{TileDB sparse arrays} in
#' \href{https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)}{COO format}:
#' \itemize{
#'  \item the row indices (\dQuote{\code{i}}) are written out as
#'   \dQuote{\code{soma_dim_0}}
#'  \item the column indices (\dQuote{\code{j}}) are written out as
#'   \dQuote{\code{soma_dim_1}}
#'  \item the non-zero values (\dQuote{\code{x}}) are written out as
#'   \dQuote{\code{soma_data}}
#' }
#' The array type is determined by \code{type}, or
#' \code{\link[arrow:infer_type]{arrow::infer_type}(slot(x, "x"))}
#'
#' @method write_soma TsparseMatrix
#' @export
#'
write_soma.TsparseMatrix <- function(
  x,
  uri,
  soma_parent,
  type = NULL,
  transpose = FALSE,
  ...,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  absolute = FALSE
) {
  stopifnot(
    "'x' must be a general sparse matrix" = inherits(x = x, what = 'generalMatrix'),
    "'x' must not be a pattern matrix" = !inherits(x = x, what = 'nsparseMatrix'),
    "'type' must be an Arrow type" = is.null(type) ||
      (R6::is.R6(type) && inherits(x = type, what = 'DataType')),
    "'transpose' must be a single logical value" = is_scalar_logical(transpose)
  )
  # Create a proper URI
  uri <- .check_soma_uri(
    uri = uri,
    soma_parent = soma_parent,
    absolute = absolute
  )
  # Transpose the matrix
  if (isTRUE(transpose)) {
    x <- Matrix::t(x)
  }
  # Create the array
  array <- SOMASparseNDArrayCreate(
    uri = uri,
    type = type %||% arrow::infer_type(methods::slot(object = x, name = 'x')),
    shape = dim(x),
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  # Write and return
  array$write(x)
  return(array)
}

#' Add an index to a data frame
#'
#' @details The index column will be determined based on availability in
#' \code{x}'s existing names. The potential names, in order of preference, are:
#' \itemize{
#'  \item \dQuote{\code{index}}
#'  \item \dQuote{\code{_index}}
#'  \item \code{alt}
#'  \item \code{paste0(prefix, "_index")}
#'  \item \code{paste(prefix, alt, sep = "_")}
#' }
#'
#' @inheritParams SeuratObject::RandomName
#' @param x A \code{data.frame}
#' @param alt An alternate index name
#' @param prefix Prefix for alternate indexes
#'
#' @return \code{x} with the row names added as an index and an attribute named
#' \dQuote{\code{index}} naming the index column
#'
#' @keywords internal
#'
#' @noRd
#'
.df_index <- function(x, alt = 'rownames', prefix = 'tiledbsoma', ...) {
  .check_seurat_installed()
  stopifnot(
    "'x' must be a data frame" = is.data.frame(x),
    "'alt' must be a single character value" = is_scalar_character(alt),
    "'prefix' must be a single character value" = is_scalar_character(prefix)
  )
  index <- ''
  i <- 1L
  while (!nzchar(index) || index %in% names(x)) {
    index <- switch(
      EXPR = i,
      '1' = 'index',
      '2' = '_index',
      '3' = alt,
      '4' = paste0(prefix, '_index'),
      '5' = paste(prefix, alt, sep = '_'),
      SeuratObject::RandomName(length = i, ...)
    )
    i <- i + 1L
  }
  x[[index]] <- row.names(x)
  attr(x = x, which = 'index') <- index
  return(x)
}

.check_soma_uri <- function(
  uri,
  soma_parent = NULL,
  absolute = FALSE,
  overwrite = FALSE
) {
  stopifnot(
    "'uri' must be a single character value" = is_scalar_character(uri),
    "'soma_parent' must be a SOMACollection" = is.null(soma_parent) ||
      inherits(x = soma_parent, what = 'SOMACollectionBase'),
    "'absolute' must be a single logical value" = is_scalar_logical(absolute)
  )
  if (!isTRUE(absolute)) {
    if (basename(uri) != uri) {
      warning("uri", call. = FALSE, immediate. = TRUE)
      uri <- basename(uri)
    }
    uri <- file_path(soma_parent$uri %||% user_dir(), uri)
  } else if (!is_remote_uri(uri)) {
    if (isTRUE(overwrite)) {
      unlink(x = uri, recursive = TRUE, force = TRUE)
    }
    dir.create(dirname(uri), recursive = TRUE)
  }
  return(uri)
}
