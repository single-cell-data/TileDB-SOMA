#' Write a SOMA Object from an \R Object
#'
#' Convert \R objects to their appropriate SOMA counterpart
#' function and methods can be written for it to provide a high-level
#' \R \eqn{\rightarrow} SOMA interface.
#'
#' @param x An object.
#' @param uri URI for resulting SOMA object.
#' @template param-dots-method
#' @param platform_config Optional \link[tiledbsoma:PlatformConfig]{platform
#' configuration}.
#' @param tiledbsoma_ctx Optional \code{\link{SOMATileDBContext}}.
#'
#' @return The URI to the resulting \code{\link{SOMAExperiment}} generated from
#' the data contained in \code{x}.
#'
#' @section Known methods:
#' \itemize{
#'  \item \link[tiledbsoma:write_soma.Seurat]{Writing Seurat objects}.
#'  \item \link[tiledbsoma:write_soma.SummarizedExperiment]{Writing SummarizedExperiment objects}.
#'  \item \link[tiledbsoma:write_soma.SingleCellExperiment]{Writing SingleCellExperiment objects}.
#' }
#'
#' @export
#'
#' @inherit write_soma_objects examples
#'
write_soma <- function(x, uri, ..., platform_config = NULL, tiledbsoma_ctx = NULL) {
  UseMethod(generic = "write_soma", object = x)
}

#' Write R Objects to SOMA
#'
#' Various helpers to write R objects to SOMA.
#'
#' @inheritParams write_soma
#' @param soma_parent The parent \link[tiledbsoma:SOMACollection]{collection}
#' (eg. a \code{\link{SOMACollection}}, \code{\link{SOMAExperiment}}, or
#' \code{\link{SOMAMeasurement}}).
#' @param ingest_mode Ingestion mode when creating the SOMA; choose from:
#' \itemize{
#'  \item \dQuote{\code{write}}: create a new SOMA and error if it already
#'   exists.
#'  \item \dQuote{\code{resume}}: attempt to create a new SOMA; if it already
#'   exists, simply open it for writing.
#' }
#' @param relative \strong{\[Internal use only\]} Is \code{uri}
#' relative or absolute.
#'
#' @return The resulting SOMA \link[tiledbsoma:SOMASparseNDArray]{array} or
#' \link[tiledbsoma:SOMADataFrame]{data frame}, returned opened for write.
#'
#' @name write_soma_objects
#' @rdname write_soma_objects
#'
#' @keywords internal
#'
NULL

#' @rdname write_soma_objects
#'
#' @section Writing Character Arrays:
#' \link[base:character]{Characters} are written out as
#' \code{\link{SOMADataFrame}s} with one attribute: \dQuote{\code{values}};
#' additionally one bit of array-level metadata is added:
#' \itemize{
#'  \item \dQuote{\code{soma_uns_outgest_hint}} with a value of
#'   \dQuote{\code{array_1d}}.
#' }
#'
#' @method write_soma character
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' # Write a character vector to a SOMA
#' uri <- withr::local_tempfile(pattern = "character")
#' (sdf <- write_soma(letters, uri, soma_parent = NULL, relative = FALSE))
#'
#' sdf$close()
#'
write_soma.character <- function(
  x,
  uri,
  soma_parent,
  ...,
  key = NULL,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  relative = TRUE
) {
  sdf <- write_soma(
    x = data.frame(values = x),
    uri = uri,
    soma_parent = soma_parent,
    df_index = "values",
    ...,
    key = key,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    relative = relative
  )
  sdf$set_metadata(uns_hint("1d"))
  return(sdf)
}

#' @param df_index The name of the column in \code{x} with the index
#' (row names); by default, will automatically add the row names of \code{x}
#' to an attribute named \dQuote{\code{index}} to the resulting
#' \code{\link{SOMADataFrame}}.
#' @param index_column_names Names of columns in \code{x} to index in the
#' resulting SOMA object.
#' @param key Optionally register the resulting \code{SOMADataFrame} in
#' \code{soma_parent} as \code{key}; pass \code{NULL} to prevent registration
#' to handle manually.
#'
#' @name write_soma_objects
#' @rdname write_soma_objects
#'
#' @section Writing Data Frames:
#' \link[base:data.frame]{Data frames} are written out as
#' \code{\link{SOMADataFrame}s}. The following transformations
#' are applied to \code{x}:
#' \itemize{
#'  \item row names are added to a column in \code{x} entitled
#'   \dQuote{\code{index}}, \dQuote{\code{_index}}, or a random name if
#'   either option is already present in \code{x}.
#'  \item a column \dQuote{\code{soma_joinid}} will be automatically
#'   added going from \code{[0, nrow(x) - 1]} encoded as
#'   \link[bit64:integer64]{64-bit integers}.
#' }
#' The array type for each column will be determined by
#' \code{\link[arrow:infer_type]{arrow::infer_type}()}; if any column contains
#' a \link[base:is.atomic]{non-atomic} type (excluding
#' \link[base:factor]{factors}, \code{\link[base]{complex}es},and
#' \code{\link[base]{raw}s}), the code will error out.
#'
#' @method write_soma data.frame
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE) && requireNamespace("SeuratObject", quietly = TRUE)
#' # Write a data.frame to a SOMA
#' uri <- withr::local_tempfile(pattern = "data-frame")
#' data("pbmc_small", package = "SeuratObject")
#' head(obs <- suppressWarnings(SeuratObject::UpdateSeuratObject(pbmc_small))[[]])
#'
#' (sdf <- write_soma(obs, uri, soma_parent = NULL, relative = FALSE))
#'
#' sdf$close()
#'
write_soma.data.frame <- function(
  x,
  uri,
  soma_parent,
  df_index = NULL,
  index_column_names = "soma_joinid",
  ...,
  key = NULL,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  relative = TRUE
) {
  stopifnot(
    "'x' must be named" = is_named(x, allow_empty = FALSE),
    "'x' must have at lease one row and one column" = dim(x) > 0L,
    "'df_index' must be a single character value" = is.null(df_index) ||
      (is_scalar_character(df_index) && nzchar(df_index)),
    "'index_column_names' must be a character vector" = is.character(index_column_names),
    "'key' must be a single character value" = is.null(key) ||
      (is_scalar_character(key) && nzchar(key))
  )
  ingest_mode <- match.arg(arg = ingest_mode, choices = c("write", "resume"))
  # Create a proper URI
  uri <- .check_soma_uri(
    uri = uri,
    soma_parent = soma_parent,
    relative = relative
  )
  if (is.character(key) && is.null(soma_parent)) {
    stop("'soma_parent' must be a SOMACollection if 'key' is provided")
  }
  # Clean up data types in `x`
  remove <- vector(mode = "logical", length = ncol(x))
  for (i in seq_len(ncol(x))) {
    col <- names(x)[i]
    remove[i] <- !inherits(
      x = try(expr = arrow::infer_type(x[[col]]), silent = TRUE),
      what = "DataType"
    )
  }
  if (any(remove)) {
    stop(
      paste(
        strwrap(paste(
          "The following columns contain unsupported data types:",
          string_collapse(sQuote(names(x)[remove]))
        )),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }
  # Check enumerations
  enumerations <- sapply(x, FUN = levels, simplify = FALSE, USE.NAMES = TRUE)
  if (all(vapply_lgl(enumerations, is.null))) {
    enumerations <- NULL
  }
  # Check `df_index`
  df_index <- df_index %||% attr(x = x, which = "index")
  if (is.null(df_index)) {
    x <- .df_index(x = x, ...)
    # x <- .df_index(x = x)
    df_index <- attr(x = x, which = "index")
  }
  if (!df_index %in% names(x)) {
    stop(
      "Unable to find ",
      sQuote(df_index),
      " in the provided data frame",
      call. = FALSE
    )
  }
  # Add `soma_joinid` to `x`
  if (!"soma_joinid" %in% names(x)) {
    # bit64::seq.integer64 does not support seq(from = 0, to = 0)
    x$soma_joinid <- if (nrow(x) == 1L) {
      bit64::integer64(length = 1L)
    } else {
      seq(bit64::as.integer64(0L), to = nrow(x) - 1L)
    }
  }
  # Check `index_column_names`
  index_column_names <- match.arg(
    arg = index_column_names,
    choices = names(x),
    several.ok = TRUE
  )

  # For index_column_name being soma_joinid -- this being the default
  # -- set that domain slot to match the data. This will endow the
  # dataframe with something users think of as a "shape". For the
  # other slots, set the domain wide open.
  #
  domain <- list()
  for (index_column_name in index_column_names) {
    if (index_column_name == "soma_joinid") {
      domain[["soma_joinid"]] <- c(0, nrow(x) - 1)
    } else {
      domain[[index_column_name]] <- NULL
    }
  }

  # Create the SOMADataFrame
  tbl <- arrow::arrow_table(x)
  sdf <- SOMADataFrameCreate(
    uri = uri,
    schema = tbl$schema,
    index_column_names = index_column_names,
    domain = domain,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  # Write values
  if (ingest_mode %in% c("resume")) {
    join_ids <- .read_soma_joinids(sdf)
    idx <- which(!x$soma_joinid %in% join_ids)
    tbl <- if (length(idx)) {
      tbl[idx, , drop = FALSE]
    } else {
      NULL
    }
  }
  if (ingest_mode %in% c("resume") && sdf$tiledbsoma_has_upgraded_domain()) {
    sdf$tiledbsoma_resize_soma_joinid_shape(nrow(x))
  }
  if (!is.null(tbl)) {
    sdf$write(tbl)
  }
  # Add to `soma_parent`
  if (is.character(key)) {
    mode <- sdf$mode()
    on.exit(sdf$reopen(mode), add = TRUE, after = FALSE)
    withCallingHandlers(
      expr = .register_soma_object(sdf, soma_parent, key, relative),
      existingKeyWarning = .maybe_muffle
    )
  }
  # Return
  return(sdf)
}

#' @param sparse Create a \link[tiledbsoma:SOMASparseNDArray]{sparse} or
#' \link[tiledbsoma:SOMADenseNDArray]{dense} array from \code{x}.
#' @param type \link[arrow:data-type]{Arrow type} for encoding \code{x}
#' (eg. \code{\link[arrow:data-type]{arrow::int32}()}); by default, attempts to
#' determine arrow type with
#' \code{\link[arrow:infer_type]{arrow::infer_type}()}.
#' @param transpose Transpose \code{x} before writing.
#' @param shape A vector of two positive integers giving the on-disk shape of
#' the array; defaults to \code{dim(x)}.
#'
#'
#' @name write_soma_objects
#' @rdname write_soma_objects
#'
#' @section Writing Dense Matrices:
#' Dense matrices are written as two-dimensional
#' \link[tiledbsoma:SOMADenseNDArray]{dense arrays}. The overall shape of the
#' array is determined by \code{dim(x)} and the type of the array is determined
#' by \code{type} or \code{\link[arrow:infer_type]{arrow::infer_type}(x)}.
#'
#' @method write_soma matrix
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' # Write a matrix to a SOMA
#' uri <- withr::local_tempfile(pattern = "matrix")
#' mat <- matrix(stats::rnorm(25L), nrow = 5L, ncol = 5L)
#' (arr <- write_soma(mat, uri, soma_parent = NULL, sparse = FALSE, relative = FALSE))
#'
#' arr$close()
#'
write_soma.matrix <- function(
  x,
  uri,
  soma_parent,
  sparse = TRUE,
  type = NULL,
  transpose = FALSE,
  ...,
  key = NULL,
  ingest_mode = "write",
  shape = NULL,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  relative = TRUE
) {
  stopifnot(
    "'sparse' must be a single logical value" = is_scalar_logical(sparse),
    "'type' must be an Arrow type" = is.null(type) || is_arrow_data_type(type),
    "'transpose' must be a single logical value" = is_scalar_logical(transpose),
    "'key' must be a single character value" = is.null(key) ||
      (is_scalar_character(key) && nzchar(key)),
    "'shape' must be a vector of two postiive integers" = is.null(shape) ||
      (rlang::is_integerish(shape, n = 2L, finite = TRUE) && all(shape > 0L))
  )
  ingest_mode <- match.arg(arg = ingest_mode, choices = c("write", "resume"))
  if (!isTRUE(sparse) && inherits(x = x, what = "sparseMatrix")) {
    stop(
      "A sparse matrix was provided and a dense array was asked for",
      call. = FALSE
    )
  }
  # Create a sparse array
  if (isTRUE(sparse)) {
    return(write_soma(
      x = methods::as(object = x, Class = "TsparseMatrix"),
      uri = uri,
      soma_parent = soma_parent,
      type = type,
      transpose = transpose,
      ...,
      key = key,
      ingest_mode = ingest_mode,
      shape = shape,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      relative = relative
    ))
  }
  if (ingest_mode != "write") {
    stop(
      "Ingestion mode of ",
      sQuote(ingest_mode, q = FALSE),
      " is not supported with dense arrays",
      call. = FALSE
    )
  }
  # Create a dense array
  if (inherits(x = x, what = "Matrix")) {
    x <- as.matrix(x)
  }
  # Create a proper URI
  uri <- .check_soma_uri(
    uri = uri,
    soma_parent = soma_parent,
    relative = relative
  )
  if (is.character(key) && is.null(soma_parent)) {
    stop("'soma_parent' must be a SOMACollection if 'key' is provided")
  }
  # Transpose the matrix
  if (isTRUE(transpose)) {
    x <- t(x)
  }
  # Check the shape
  if (!is.null(shape) && any(shape < dim(x))) {
    stop(
      "Requested an array of shape (",
      paste(shape, collapse = ", "),
      "), but was given a matrix with a larger shape (",
      paste(dim(x), collapse = ", "),
      ")",
      call. = FALSE
    )
  }
  # Create the array
  array <- SOMADenseNDArrayCreate(
    uri = uri,
    type = type %||% arrow::infer_type(x),
    shape = shape %||% dim(x),
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  # Write values
  array$write(x)
  # Add to `soma_parent`
  if (is.character(key)) {
    mode <- array$mode()
    on.exit(array$reopen(mode), add = TRUE, after = FALSE)
    withCallingHandlers(
      expr = .register_soma_object(array, soma_parent, key, relative),
      existingKeyWarning = .maybe_muffle
    )
  }
  # Return
  return(array)
}

#' @name write_soma_objects
#' @rdname write_soma_objects
#'
#' @method write_soma Matrix
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' # Write a dense S4 Matrix to a SOMA
#' uri <- withr::local_tempfile(pattern = "s4-matrix")
#' mat <- Matrix::Matrix(stats::rnorm(25L), nrow = 5L, ncol = 5L)
#' (arr <- write_soma(mat, uri, soma_parent = NULL, sparse = FALSE, relative = FALSE))
#'
#' arr$close()
#'
write_soma.Matrix <- write_soma.matrix

#' @name write_soma_objects
#' @rdname write_soma_objects
#'
#' @section Writing Sparse Matrices:
#' Sparse matrices are written out as two-dimensional
#' \link[tiledbsoma:SOMASparseNDArray]{TileDB sparse arrays} in
#' \href{https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)}{COO format}:
#' \itemize{
#'  \item the row indices (\dQuote{\code{i}}) are written out as
#'   \dQuote{\code{soma_dim_0}}.
#'  \item the column indices (\dQuote{\code{j}}) are written out as
#'   \dQuote{\code{soma_dim_1}}.
#'  \item the non-zero values (\dQuote{\code{x}}) are written out as
#'   \dQuote{\code{soma_data}}.
#' }
#' The array type is determined by \code{type}, or
#' \code{\link[arrow:infer_type]{arrow::infer_type}(slot(x, "x"))}.
#'
#' @method write_soma TsparseMatrix
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' # Write a TsparseMatrix to a SOMA
#' uri <- withr::local_tempfile(pattern = "tsparse-matrix")
#' mat <- Matrix::rsparsematrix(5L, 5L, 0.3, repr = "T")
#' (arr <- write_soma(mat, uri, soma_parent = NULL, relative = FALSE))
#'
#' arr$close()
#'
#' # Write a CsparseMatrix to a SOMA
#' uri <- withr::local_tempfile(pattern = "csparse-matrix")
#' mat <- Matrix::rsparsematrix(5L, 5L, 0.3, repr = "C")
#' (arr <- write_soma(mat, uri, soma_parent = NULL, relative = FALSE))
#'
#' arr$close()
#'
#' # Write an RsparseMatrix to a SOMA
#' uri <- withr::local_tempfile(pattern = "rsparse-matrix")
#' mat <- Matrix::rsparsematrix(5L, 5L, 0.3, repr = "R")
#' (arr <- write_soma(mat, uri, soma_parent = NULL, relative = FALSE))
#'
#' arr$close()
#'
write_soma.TsparseMatrix <- function(
  x,
  uri,
  soma_parent,
  type = NULL,
  transpose = FALSE,
  ...,
  key = NULL,
  ingest_mode = "write",
  shape = NULL,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  relative = TRUE
) {
  stopifnot(
    "'x' must be a general sparse matrix" = inherits(x = x, what = "generalMatrix"),
    "'x' must not be a pattern matrix" = !inherits(x = x, what = "nsparseMatrix"),
    "'type' must be an Arrow type" = is.null(type) ||
      (R6::is.R6(type) && inherits(x = type, what = "DataType")),
    "'transpose' must be a single logical value" = is_scalar_logical(transpose),
    "'key' must be a single character value" = is.null(key) ||
      (is_scalar_character(key) && nzchar(key)),
    "'shape' must be a vector of two postiive integers" = is.null(shape) ||
      (rlang::is_integerish(shape, n = 2L, finite = TRUE) && all(shape > 0L))
  )
  ingest_mode <- match.arg(arg = ingest_mode, choices = c("write", "resume"))
  # Create a proper URI
  uri <- .check_soma_uri(
    uri = uri,
    soma_parent = soma_parent,
    relative = relative
  )
  if (is.character(key) && is.null(soma_parent)) {
    stop("'soma_parent' must be a SOMACollection if 'key' is provided")
  }
  # Transpose the matrix
  if (isTRUE(transpose)) {
    x <- Matrix::t(x)
  }
  # Check the shape
  if (!is.null(shape) && any(shape < dim(x))) {
    stop(
      "Requested an array of shape (",
      paste(shape, collapse = ", "),
      "), but was given a matrix with a larger shape (",
      paste(dim(x), collapse = ", "),
      ")",
      call. = FALSE
    )
  }
  # Create the array
  array <- SOMASparseNDArrayCreate(
    uri = uri,
    type = type %||% arrow::infer_type(methods::slot(object = x, name = "x")),
    shape = shape %||% dim(x),
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    tiledb_timestamp = Sys.time()
  )
  # Write values
  if (ingest_mode %in% c("resume")) {
    if (array$ndim() != 2L) {
      stop(
        "Attempting to resume writing a matrix to a sparse array with more than two dimensions",
        call. = FALSE
      )
    }
    row_ids <- .read_soma_joinids(array, axis = 0L)
    col_ids <- .read_soma_joinids(array, axis = 1L)
    tbl <- data.frame(
      i = bit64::as.integer64(methods::slot(x, "i")),
      j = bit64::as.integer64(methods::slot(x, "j")),
      x = methods::slot(x, "x")
    )
    tbl <- tbl[-which(tbl$i %in% row_ids & tbl$j %in% col_ids), , drop = FALSE]
    x <- if (nrow(tbl)) {
      Matrix::sparseMatrix(
        i = as.integer(tbl$i),
        j = as.integer(tbl$j),
        x = tbl$x,
        dims = dim(x),
        index1 = FALSE,
        repr = "T"
      )
    } else {
      NULL
    }
  }
  if (!is.null(x)) {
    array$write(x)
  }
  # Add to `soma_parent`
  if (is.character(key)) {
    mode <- array$mode()
    on.exit(array$reopen(mode), add = TRUE, after = FALSE)
    withCallingHandlers(
      expr = .register_soma_object(array, soma_parent, key, relative),
      existingKeyWarning = .maybe_muffle
    )
  }
  # Return
  return(array)
}

#' Add an index to a data frame
#'
#' @details The index column will be determined based on availability in
#' \code{x}'s existing names. The potential names, in order of preference, are:
#' \itemize{
#'  \item \dQuote{\code{obs_id}} or \dQuote{\code{var_id}},
#'    depending on \code{axis}.
#'  \item \code{alt}.
#'  \item \code{paste0(prefix, "obs_id")} or \code{paste0(prefix, "var_id")}.
#'  \item \code{paste(prefix, alt, sep = "_")}.
#' }
#'
#' @inheritParams SeuratObject::RandomName
#' @param x A \code{data.frame}.
#' @param alt An alternate index name.
#' @param axis Either \dQuote{\code{obs}} or \dQuote{\code{var}} for
#' default index name.
#' @param prefix Prefix for alternate indexes.
#'
#' @return \code{x} with the row names added as an index and an attribute named
#' \dQuote{\code{index}} naming the index column.
#'
#' @keywords internal
#'
#' @noRd
#'
.df_index <- function(
  x,
  alt = "rownames",
  axis = "obs",
  prefix = "tiledbsoma",
  ...
) {
  stopifnot(
    "'x' must be a data frame" = is.data.frame(x) || inherits(x, "DataFrame"),
    "'alt' must be a single character value" = is_scalar_character(alt),
    "'axis' must be a single character value" = is_scalar_character(axis),
    "'prefix' must be a single character value" = is_scalar_character(prefix)
  )
  axis <- match.arg(axis, choices = c("obs", "var", "index"))
  default <- switch(EXPR = axis,
    index = "index",
    paste0(axis, "_id")
  )
  index <- ""
  i <- 1L
  while (!nzchar(index) || index %in% names(x)) {
    index <- switch(
      EXPR = i,
      "1" = default,
      "2" = alt,
      "3" = paste(prefix, default, sep = "_"),
      "4" = paste(prefix, alt, sep = "_"),
      random_name(length = i, ...)
    )
    i <- i + 1L
  }
  x[[index]] <- row.names(x)
  attr(x = x, which = "index") <- index
  return(x)
}

.check_soma_uri <- function(
  uri,
  soma_parent = NULL,
  relative = TRUE
) {
  stopifnot(
    "'uri' must be a single character value" = is_scalar_character(uri),
    "'soma_parent' must be a SOMACollection" = is.null(soma_parent) ||
      inherits(x = soma_parent, what = "SOMACollectionBase"),
    "'relative' must be a single logical value" = is_scalar_logical(relative)
  )
  if (!isFALSE(relative)) {
    if (basename(uri) != uri) {
      warning("uri", call. = FALSE, immediate. = TRUE)
      uri <- basename(uri)
    }
    uri <- file_path(soma_parent$uri %||% tools::R_user_dir("tiledbsoma"), uri)
  } else if (!is_remote_uri(uri)) {
    dir.create(dirname(uri), showWarnings = FALSE, recursive = TRUE)
  }
  return(uri)
}

.register_soma_object <- function(x, soma_parent, key, relative = TRUE) {
  stopifnot(
    "'x' must be a SOMA object" = inherits(x, c("SOMAArrayBase", "SOMACollectionBase")),
    "'soma_parent' must be a SOMA collection" = inherits(soma_parent, "SOMACollectionBase"),
    "'key' must be a single character value" = is_scalar_character(key) && nzchar(key),
    "'relative' must be a single logical value" = is_scalar_logical(relative)
  )
  xmode <- x$mode()
  if (xmode == "CLOSED") {
    x$reopen("READ", tiledb_timestamp = x$tiledb_timestamp)
    xmode <- x$mode()
  }
  on.exit(x$reopen(mode = xmode), add = TRUE, after = FALSE)
  oldmode <- soma_parent$mode()
  if (oldmode == "CLOSED") {
    soma_parent$reopen("READ", tiledb_timestamp = soma_parent$tiledb_timestamp)
    oldmode <- soma_parent$mode()
  }
  on.exit(soma_parent$reopen(oldmode), add = TRUE, after = FALSE)
  if (key %in% soma_parent$names()) {
    existing <- soma_parent$get(key)
    warning(warningCondition(
      message = paste(
        "Already found a",
        existing$class(),
        "stored as",
        sQuote(key),
        "in the parent collection"
      ),
      class = "existingKeyWarning"
    ))
    return(invisible(NULL))
  }
  soma_parent$reopen("WRITE")
  soma_parent$set(
    x,
    name = key,
    relative = switch(uri_scheme(x$uri) %||% "",
      tiledb = FALSE,
      relative
    )
  )
  return(invisible(NULL))
}
