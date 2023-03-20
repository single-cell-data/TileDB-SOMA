#' @importFrom rlang is_integerish
#'
NULL

#' Write a SOMA from an R Object
#'
#' @param x An object
#' @param uri URI for resulting SOMA
#' @template param-dots-method
#' @param platform_config ...
#' @param tiledbsoma_ctx ...
#'
#' @return A \code{\link{SOMAExperiment}} with the data from \code{x}
#'
#' @export
#'
write_soma <- function(x, uri, ..., platform_config = NULL, tiledbsoma_ctx = NULL) {
  UseMethod(generic = 'write_soma', object = x)
}

#' @method write_soma Assay
#' @export
#'
write_soma.Assay <- function(
  x,
  uri = NULL,
  soma,
  ...,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  absolute = FALSE
) {
  .check_seurat_installed()
  stopifnot(
    "'uri' must be a single character value" = is.null(uri) ||
      is_scalar_character(uri),
    "'soma' must be a SOMACollection" = inherits(soma, what = 'SOMACollectionBase'),
    "'absolute' must be a single logical value" = is_scalar_logical(absolute)
  )
  # Create a proper URI
  uri <- uri %||% gsub(pattern = '_$', replacement = '', x = SeuratObject::Key(x))
  uri <- .check_soma_uri(uri = uri, soma = soma, absolute = absolute)
  # Create the measurement
  ms <- SOMAMeasurementCreate(
    uri = uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  ms$X <- SOMACollectionCreate(
    uri = file_path(ms$uri, 'X'),
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  # Write `X` matrices
  for (slot in c('counts', 'data', 'scale.data')) {
    mat <- SeuratObject::GetAssayData(object = x, slot = slot)
    if (SeuratObject::IsMatrixEmpty(mat)) {
      next
    }
    if (!identical(x = dim(mat), y = dim(x))) {
      mat <- pad_matrix(
        x = SeuratObject::as.sparse(x = mat),
        rownames = rownames(x),
        colnames = colnames(x)
      )
    }
    lyr <- gsub(pattern = '\\.', replacement = '_', x = slot)
    tryCatch(
      expr = ms$X$set(
        object = write_soma(
          x = mat,
          uri = lyr,
          soma = ms$X,
          sparse = TRUE,
          transpose = TRUE,
          platform_config = platform_config,
          tiledbsoma_ctx = tiledbsoma_ctx
        ),
        name = lyr
      ),
      error = function(err) {
        if (slot == 'data') {
          stop(err)
        }
        err_to_warn(err)
      }
    )
  }
  # Write feature-level meta data
  meta_data <- .df_index(x = x[[]], alt = 'features', prefix = 'seurat')
  meta_data[[attr(x = meta_data, which = 'index')]] <- rownames(x)
  ms$var <- write_soma(
    x = meta_data,
    uri = 'var',
    soma = ms,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  # Return
  if (class(x)[1L] != 'Assay') {
    warning(
      paste(
        strwrap(paste0(
          "Extended assays (eg. ",
          class(x)[1L],
          ") cannot be written to SOMAs at this time"
        )),
        collapse = '\n'
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  return(ms)
}

#' @method write_soma data.frame
#' @export
#'
write_soma.data.frame <- function(
  x,
  uri,
  soma,
  index = NULL,
  ...,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  absolute = FALSE
) {
  stopifnot(
    "'index' must be a single character value" = is.null(index) ||
      is_scalar_character(index)
  )
  # Create a proper URI
  uri <- .check_soma_uri(uri = uri, soma = soma, absolute = absolute)
  # Clean up data types in `x`
  remove <- vector(mode = 'logical', length = ncol(x))
  for (i in seq_len(ncol(x))) {
    col <- names(x)[i]
    if (is.factor(x[[col]])) {
      x[[col]] <- as.character(x[[col]])
    }
    remove[i] <- !inherits(
      x = x[[col]],
      what = setdiff(x = .SCALAR_TYPES(), y = 'any')
    )
  }
  if (any(remove)) {
    warning("remove")
    x <- x[, !remove, drop = FALSE]
  }
  # Check `index`
  index <- index %||% attr(x = x, which = 'index')
  if (is.null(index)) {
    x <- .df_index(x = x, ...)
    index <- attr(x = x, which = 'index')
  }
  if (!index %in% names(x)) {
    stop(
      "Unable to find ",
      sQuote(index),
      " in the provided data frame",
      call. = FALSE
    )
  }
  # Add `soma_joinid` to `x`
  if ('soma_joinid' %in% names(x)) {
    warning('soma_joinid')
  }
  x$soma_joinid <- bit64::seq.integer64(from = 0L, to = nrow(x) - 1L)
  # Create the SOMADataFrame
  tbl <- arrow::arrow_table(x)
  sdf <- SOMADataFrameCreate(
    uri = uri,
    schema = tbl$schema,
    index_column_names = index,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  # Write and return
  sdf$write(tbl)
  return(sdf)
}

#' @method write_soma DimReduc
#' @export
#'
write_soma.DimReduc <- function(
  x,
  uri = NULL,
  soma,
  fidx = NULL,
  nfeatures = NULL,
  ...,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  absolute = FALSE
) {
  .check_seurat_installed()
  stopifnot(
    "'uri' must be NULL" = is.null(uri),
    "'soma' must be a SOMAMeasurement" = inherits(soma, what = 'SOMAMeasurement'),
    "'fidx' must be a positive integer vector" = is.null(fidx) ||
      (rlang::is_integerish(fidx, finite = TRUE) && all(fidx > 0L)),
    "'nfeatures' must be a single positive integer" = is.null(nfeatures) ||
      (rlang::is_integerish(nfeatures, n = 1L, finite = TRUE) && nfeatures > 0L),
    "'absolute' must be a single logical value" = is_scalar_logical(absolute)
  )
  key <- tolower(gsub(pattern = '_$', replacement = '', x = SeuratObject::Key(x)))
  key <- switch(EXPR = key, pc = 'pca', ic = 'ica', key)
  # Create a group for `obs,`
  if (!'obsm' %in% soma$names()) {
    soma$obsm <- SOMACollectionCreate(
      uri = file_path(soma$uri, 'obsm'),
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  }
  embed <- paste0('X_', key)
  soma$obsm$set(
    object = write_soma(
      x = SeuratObject::Embeddings(x),
      uri = embed,
      soma = soma$obsm,
      sparse = FALSE,
      transpose = FALSE,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    ),
    name = embed
  )
  # Add feature loadings
  loadings <- SeuratObject::Loadings(x)
  # Check feature info
  if (!SeuratObject::IsMatrixEmpty(loadings)) {
    finfo <- vapply_lgl(X = list(fidx, nfeatures), FUN = is.null)
    msg <- if (all(finfo)) {
      "No feature information provided, not adding feature loadings"
    } else if (any(finfo) && !all(finfo)) {
      paste(
        "Either both",
        sQuote('fidx'),
        "and",
        sQuote('nfeatures'),
        "must be supplied or both must be NULL"
      )
    } else if (max(fidx) > nfeatures) {
      paste(sQuote('fidx'), 'exceeds', sQuote('nfeatures'))
    } else if (all(is.na(fidx))) {
      "No feature index match"
    } else {
      ''
    }
    if (nzchar(msg)) {
      warning(
        paste(
          strwrap(paste0(msg, ', not adding feature loadings')),
          collapse = '\n'
        ),
        call. = FALSE,
        immediate. = TRUE
      )
      loadings <- methods::new('matrix')
    }
  }
  # Write feature loadings
  if (!SeuratObject::IsMatrixEmpty(loadings)) {
    ldgs <- switch(EXPR = key, pca = 'PCs', ica = 'ICs', paste0(toupper(key), 's'))
    # Create a group for `varm`
    if (!'varm' %in% soma$names()) {
      soma$varm <- SOMACollectionCreate(
        uri = file_path(soma$uri, 'varm'),
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      )
    }
    # Pad our feature loadings matrix
    mat <- matrix(data = NA_real_, nrow = nfeatures, ncol = ncol(loadings))
    mat[fidx, ] <- loadings
    # Write the feature loadings
    soma$varm$set(
      object = write_soma(
        x = mat,
        uri = ldgs,
        soma = soma$varm,
        sparse = FALSE,
        transpose = FALSE,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      ),
      name = ldgs
    )
  }
  return(invisible(soma))
}

#' @method write_soma Graph
#' @export
#'
write_soma.Graph <- function(
  x,
  uri = NULL,
  soma,
  ...,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  absolute = FALSE
) {
  .check_seurat_installed()
  stopifnot(
    "'uri' must be a single character value" = is.null(uri) ||
      is_scalar_character(uri),
    "'soma' must be a SOMAMeasurement" = inherits(soma, what = 'SOMAMeasurement'),
    "'absolute' must be a single logical value" = is_scalar_logical(absolute)
  )
  if (!'obsp' %in% soma$names()) {
    soma$obsp <- SOMACollectionCreate(
      uri = file_path(soma$uri, 'obsp'),
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  }
  soma$obsp$set(
    object = NextMethod(
      generic = 'write_soma',
      object = x,
      uri = uri,
      soma = soma$obsp,
      sparse = TRUE,
      transpose = FALSE,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    ),
    name = uri
  )
  return(invisible(soma))
}

#' @method write_soma matrix
#' @export
#'
write_soma.matrix <- function(
  x,
  uri,
  soma,
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
      soma = soma,
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
  uri <- .check_soma_uri(uri = uri, soma = soma, absolute = absolute)
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

#' @method write_soma Matrix
#' @export
#'
write_soma.Matrix <- write_soma.matrix

#' @method write_soma TsparseMatrix
#' @export
#'
write_soma.TsparseMatrix <- function(
  x,
  uri,
  soma,
  type = NULL,
  transpose = FALSE,
  ...,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  absolute = FALSE
) {
  stopifnot(
    "'type' must be an Arrow type" = is.null(type) ||
      (R6::is.R6(type) && inherits(x = type, what = 'DataType')),
    "'transpose' must be a single logical value" = is_scalar_logical(transpose)
  )
  # Create a proper URI
  uri <- .check_soma_uri(uri = uri, soma = soma, absolute = absolute)
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

#' Write a \code{\link[SeuratObject]{Seurat}} object to a SOMA
#'
#' @inheritParams write_soma
#' @param x A \code{\link[SeuratObject]{Seurat}} object
#'
#' @inherit write_soma return
#'
#' @method write_soma Seurat
#' @export
#'
write_soma.Seurat <- function(
  x,
  uri = NULL,
  ...,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  overwrite = FALSE
) {
  .check_seurat_installed()
  stopifnot(
    "'uri' must be a single character value" = is.null(uri) ||
      is_scalar_character(uri)
  )
  uri <- uri %||% file_path(user_dir(), SeuratObject::Project(x))
  if (!is_remote_uri(uri)) {
    if (isTRUE(overwrite)) {
      unlink(x = uri, recursive = TRUE, force = TRUE)
    }
    dir.create(dirname(uri), recursive = TRUE)
  }
  experiment <- SOMAExperimentCreate(
    uri = uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  # Write cell-level meta data
  meta_data <- .df_index(x = x[[]], alt = 'cells', prefix = 'seurat')
  meta_data[[attr(meta_data, 'index')]] <- colnames(x)
  experiment$obs <- write_soma(
    x = meta_data,
    uri = 'obs',
    soma = experiment,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  # Write assays
  experiment$add_new_collection(
    object = SOMACollectionCreate(
      uri = file_path(experiment$uri, 'ms'),
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    ),
    key = 'ms'
  )
  for (assay in SeuratObject::Assays(x)) {
    tryCatch(
      expr = experiment$ms$set(
        object = write_soma(
          x = x[[assay]],
          uri = assay,
          soma = experiment$ms,
          platform_config = platform_config,
          tiledbsoma_ctx = tiledbsoma_ctx
        ),
        name = assay
      ),
      error = function(err) {
        if (assay == SeuratObject::DefaultAssay(x)) {
          stop(err)
        }
        err_to_warn(err)
      }
    )
  }
  # Write dimensional reductions
  for (reduc in SeuratObject::Reductions(x)) {
    assay <- SeuratObject::DefaultAssay(x[[reduc]])
    ms <- if (assay %in% experiment$ms$names()) {
      experiment$ms$get(assay)
    } else if (SeuratObject::IsGlobal(x[[reduc]])) {
      assay <- SeuratObject::DefaultAssay(x)
      warning(
        paste(
          strwrap(paste0(
            "Cannot find a measurement for global reduction ",
            sQuote(reduc),
            " (default assay: ",
            sQuote(SeuratObject::DefaultAssay(x[[reduc]])),
            "), adding to measurement for ",
            sQuote(assay)
          )),
          collapse = '\n'
        ),
        call. = FALSE,
        immediate. = TRUE
      )
      experiment$ms$get(assay)
    } else {
      # This should never happen
      warning(
        paste(
          strwrap(paste0(
            "Cannot find a measurement for non-global reduction ",
            sQuote(reduc),
            " (default assay: ",
            sQuote(assay),
            "), skipping"
          )),
          collapse = '\n'
        ),
        call. = FALSE,
        immediate. = TRUE
      )
      next
    }
    loadings <- SeuratObject::Loadings(x[[reduc]])
    if (!SeuratObject::IsMatrixEmpty(loadings)) {
      fidx <- match(x = rownames(loadings), table = rownames(x[[assay]]))
      nfeatures <- nrow(x[[assay]])
    } else {
      fidx <- nfeatures <- NULL
    }
    tryCatch(
      expr = write_soma(
        x = x[[reduc]],
        uri = NULL,
        soma = ms,
        fidx = fidx,
        nfeatures = nfeatures,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      ),
      error = err_to_warn
    )
  }
  # Write graphs
  for (graph in SeuratObject::Graphs(x)) {
    assay <- SeuratObject::DefaultAssay(x[[graph]])
    if (!assay %in% experiment$ms$names()) {
      warning(
        paste(
          strwrap(paste0(
            "Cannot find a measurement for graph ",
            sQuote(graph),
            " (default assay: ",
            sQuote(assay),
            "), skipping"
          )),
          collapse = FALSE
        ),
        call. = FALSE,
        immediate. = TRUE
      )
      next
    }
    tryCatch(
      expr = write_soma(
        x = x[[graph]],
        uri = graph,
        soma = experiment$ms$get(assay),
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      ),
      error = err_to_warn
    )
  }
  # TODO: Write images
  if (length(SeuratObject::Images(x))) {
    warning(
      "Spatially resolved data cannot be written to SOMAs at this time",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  return(experiment)
}

#' Add an index to a data frame
#'
#' @details The index column will be determined based on availability in
#' \code{x}'s existing names. The potential names, in order of preference, are:
#' \itemize{
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
.df_index <- function(x, alt, prefix = 'tiledbsoma', ...) {
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
      '1' = '_index',
      '2' = alt,
      '3' = paste0(prefix, '_index'),
      '4' = paste(prefix, alt, sep = '_'),
      SeuratObject::RandomName(length = i, ...)
    )
    i <- i + 1L
  }
  x[[index]] <- row.names(x)
  attr(x = x, which = 'index') <- index
  return(x)
}

.check_soma_uri <- function(uri, soma = NULL, absolute = FALSE, overwrite = FALSE) {
  stopifnot(
    "'uri' must be a single character value" = is_scalar_character(uri),
    "'soma' must be a SOMACollection" = is.null(soma) ||
      inherits(x = soma, what = 'SOMACollectionBase'),
    "'absolute' must be a single logical value" = is_scalar_logical(absolute)
  )
  if (!isTRUE(absolute)) {
    if (basename(uri) != uri) {
      warning("uri", call. = FALSE, immediate. = TRUE)
      uri <- basename(uri)
    }
    uri <- file_path(soma$uri %||% user_dir(), uri)
  } else if (!is_remote_uri(uri)) {
    if (isTRUE(overwrite)) {
      unlink(x = uri, recursive = TRUE, force = TRUE)
    }
    dir.create(dirname(uri), recursive = TRUE)
  }
  return(uri)
}
