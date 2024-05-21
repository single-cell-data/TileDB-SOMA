#' @importFrom rlang is_integerish
#'
NULL

#' Convert a \pkg{Seurat} Sub-Object to a SOMA Object, returned opened for write
#'
#' Various helpers to write \pkg{Seurat} sub-objects to SOMA objects.
#'
#' @inheritParams write_soma_objects
#' @param x A \pkg{Seurat} sub-object (eg. an
#' \code{\link[SeuratObject]{Assay}}, a \code{\link[SeuratObject]{DimReduc}},
#' or a \code{\link[SeuratObject]{Graph}})
#' @param soma_parent The parent \link[tiledbsoma:SOMACollection]{collection};
#' for the \code{DimReduc} and \code{Graph} methods, this \strong{must} be a
#' \link[tiledbsoma:SOMAMeasurement]{measurement} for the assay \code{x}
#' was generated from
#'
#' @name write_soma_seurat_sub_objects
#' @rdname write_soma_seurat_sub
#'
#' @keywords internal
#'
NULL

#' @return \code{Assay} method: a \code{\link{SOMAMeasurement}} with the
#' data from \code{x}, returned opened for write
#'
#' @rdname write_soma_seurat_sub
#'
#' @section Writing \code{\link[SeuratObject]{Assay}s}:
#' \pkg{Seurat} \code{\link[SeuratObject]{Assay}} objects are written out as
#' individual \link[tiledbsoma:SOMAMeasurement]{measurements}:
#' \itemize{
#'  \item the \dQuote{\code{data}} matrix is written out a
#'   \link[tiledbsoma:SOMASparseNDArray]{sparse matrix} called
#'   \dQuote{\code{data}} within the \dQuote{\code{X}} group
#'  \item the \dQuote{\code{counts}} matrix, if not
#'   \link[SeuratObject:IsMatrixEmpty]{empty}, is written out a
#'   \link[tiledbsoma:SOMASparseNDArray]{sparse matrix} called
#'   \dQuote{\code{counts}} within the \dQuote{\code{X}} group
#'  \item the \dQuote{\code{scale.data}} matrix, if not
#'   \link[SeuratObject:IsMatrixEmpty]{empty}, is written out a
#'   \link[tiledbsoma:SOMASparseNDArray]{sparse matrix} called
#'   \dQuote{\code{scale_data}} within the \dQuote{\code{X}} group
#'  \item feature-level meta data is written out as a
#'   \link[tiledbsoma:SOMADataFrame]{data frame} called \dQuote{\code{var}}
#' }
#' Expression matrices are transposed (cells as rows) prior to writing. All
#' other slots, including results from extended assays (eg. \code{SCTAssay},
#' \code{ChromatinAssay}) are lost
#'
#' @method write_soma Assay
#' @export
#'
write_soma.Assay <- function(
  x,
  uri = NULL,
  soma_parent,
  ...,
  ingest_mode = 'write',
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  relative = TRUE
) {
  check_package('SeuratObject', version = .MINIMUM_SEURAT_VERSION())
  stopifnot(
    "'uri' must be a single character value" = is.null(uri) ||
      is_scalar_character(uri),
    "'soma_parent' must be a SOMACollection" = inherits(
      x = soma_parent,
      what = 'SOMACollectionBase'
    ),
    "'relative' must be a single logical value" = is_scalar_logical(relative)
  )

  # Find `shape` if and only if we're called from `write_soma.Seurat()`
  parents <- unique(sys.parents())
  idx <- which(vapply_lgl(
    parents,
    FUN = function(i) identical(sys.function(i), write_soma.Seurat)
  ))
  shape <- if (length(idx) == 1L) {
    get("shape", envir = sys.frame(parents[idx]))
  } else {
    NULL
  }
  shape <- rev(shape)

  # Create a proper URI
  uri <- uri %||% gsub(pattern = '_$', replacement = '', x = SeuratObject::Key(x))
  uri <- .check_soma_uri(
    uri = uri,
    soma_parent = soma_parent,
    relative = relative
  )

  # Create the measurement
  ms <- SOMAMeasurementCreate(
    uri = uri,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  X <- if (!'X' %in% ms$names()) {
    SOMACollectionCreate(
      uri = file_path(ms$uri, 'X'),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  } else if (isTRUE(relative)) {
    SOMACollectionOpen(uri = file_path(ms$uri, 'X'), mode = 'WRITE')
  } else {
    ms$X
  }
  withCallingHandlers(
    .register_soma_object(X, soma_parent = ms, key = 'X', relative = relative),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(X$close(), add = TRUE, after = FALSE)

  # Write `X` matrices
  for (slot in c("counts", "data", "scale.data")) {
    mat <- SeuratObject::GetAssayData(x, slot)
    if (SeuratObject::IsMatrixEmpty(mat)) next

    # Skip 'data' slot if it's identical to 'counts'
    if (slot == "data") {
      if (identical(mat, SeuratObject::GetAssayData(x, "counts"))) {
        spdl::info("Skipping 'data' slot because it's identical to 'counts'")
        next
      }
    }

    if (!identical(x = dim(mat), y = dim(x))) {
      spdl::info("Padding layer '{}' to match dimensions of assay", slot)
      mat <- pad_matrix(
        x = mat,
        rowidx = match(x = rownames(mat), table = rownames(x)),
        colidx = match(x = colnames(mat), table = colnames(x)),
        shape = dim(x),
        sparse = TRUE,
        rownames = rownames(x),
        colnames = colnames(x)
      )
    }

    layer <- gsub(pattern = '\\.', replacement = '_', x = slot)
    spdl::info("Adding '{}' matrix as '{}'", slot, layer)
    tryCatch(
      expr = write_soma(
        x = mat,
        uri = layer,
        soma_parent = X,
        sparse = TRUE,
        transpose = TRUE,
        ingest_mode = ingest_mode,
        shape = shape,
        key = layer,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
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
  var_df <- .df_index(
    x = x[[]],
    alt = 'features',
    axis = 'var',
    prefix = 'seurat'
  )
  var_df[[attr(x = var_df, which = 'index')]] <- rownames(x)
  spdl::info("Adding feature-level meta data")
  write_soma(
    x = var_df,
    uri = 'var',
    soma_parent = ms,
    key = 'var',
    ingest_mode = ingest_mode,
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
          ") are not fully supported; core Assay data has been written but ",
          "additional slots have been skipped"
        )),
        collapse = '\n'
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  return(ms)
}

#' @param fidx An integer vector describing the location of features in
#' \code{SeuratObject::Loadings(x)} with relation to \code{soma_parent}
#' (eg. \code{match(rownames(Loadings(x)), rownames(assay))})
#' @param nfeatures The number of features present in \code{soma_parent}
#' (eg. \code{nrow(assay)})
#'
#' @return \code{DimReduc} and \code{Graph} methods: invisibly returns
#' \code{soma_parent}, opened for write, with the values of \code{x} added to it
#'
#' @rdname write_soma_seurat_sub
#'
#' @section Writing \code{\link[SeuratObject]{DimReduc}s}:
#' \pkg{Seurat} \code{\link[SeuratObject]{DimReduc}} objects are written out
#' to the \dQuote{\code{obsm}} and \dQuote{\code{varm}} groups of a
#' \link[tiledbsoma:SOMAMeasurement]{measurement}:
#' \itemize{
#'  \item cell embeddings are written out as a
#'   \link[tiledbsoma:SOMASparseNDArray]{sparse matrix} in the
#'   \dQuote{\code{obsm}} group
#'  \item feature loadings, if not \link[SeuratObject:IsMatrixEmpty]{empty},
#'   are written out as a \link[tiledbsoma:SOMASparseNDArray]{sparse matrix} in
#'   the \dQuote{\code{varm}} groups; loadings are padded with \code{NAs}
#'   to include all features
#' }
#' Dimensional reduction names are translated to AnnData-style names (eg.
#' \dQuote{\code{pca}} becomes \code{X_pca} for embeddings and
#' \dQuote{\code{PCs}} for loadings). All other slots, including projected
#' feature loadings and jackstraw information, are lost
#'
#' @method write_soma DimReduc
#' @export
#'
write_soma.DimReduc <- function(
  x,
  uri = NULL,
  soma_parent,
  fidx = NULL,
  nfeatures = NULL,
  ...,
  ingest_mode = 'write',
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  relative = TRUE
) {
  check_package('SeuratObject', version = .MINIMUM_SEURAT_VERSION())
  stopifnot(
    "'uri' must be NULL" = is.null(uri),
    "'soma_parent' must be a SOMAMeasurement" = inherits(
      x = soma_parent,
      what = 'SOMAMeasurement'
    ),
    "'fidx' must be a positive integer vector" = is.null(fidx) ||
      (rlang::is_integerish(fidx, finite = TRUE) && all(fidx > 0L)),
    "'nfeatures' must be a single positive integer" = is.null(nfeatures) ||
      (rlang::is_integerish(nfeatures, n = 1L, finite = TRUE) && nfeatures > 0L),
    "'relative' must be a single logical value" = is_scalar_logical(relative)
  )

  key <- tolower(gsub(pattern = '_$', replacement = '', x = SeuratObject::Key(x)))
  key <- switch(EXPR = key, pc = 'pca', ic = 'ica', key)

  # Find `shape` if and only if we're called from `write_soma.Seurat()`
  parents <- unique(sys.parents())
  idx <- which(vapply_lgl(
    parents,
    FUN = function(i) identical(sys.function(i), write_soma.Seurat)
  ))
  shape <- if (length(idx) == 1L) {
    get("shape", envir = sys.frame(parents[idx]))
  } else {
    NULL
  }

  if (!is.null(shape)) {
    demb <- c(shape[2L], ncol(x))
    dload <- c(shape[1L], ncol(x))
  } else {
    demb <- dload <- NULL
  }

  # Create a group for `obsm,`
  obsm <- if (!'obsm' %in% soma_parent$names()) {
    SOMACollectionCreate(
      uri = file_path(soma_parent$uri, 'obsm'),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  } else if (isTRUE(relative)) {
    SOMACollectionOpen(uri = file_path(soma_parent$uri, 'obsm'), mode = 'WRITE')
  } else {
    soma_parent$obsm
  }
  withCallingHandlers(
    .register_soma_object(obsm, soma_parent, key = 'obsm', relative = relative),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(obsm$close(), add = TRUE, after = FALSE)

  embed <- paste0('X_', key)
  spdl::info("Adding embeddings as {}", sQuote(embed))

  # Always write reductions as sparse arrays
  write_soma(
    x = SeuratObject::Embeddings(x),
    uri = embed,
    soma_parent = obsm,
    sparse = TRUE,
    transpose = FALSE,
    key = embed,
    ingest_mode = ingest_mode,
    shape = demb,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )

  # Add feature loadings
  loadings <- SeuratObject::Loadings(x)

  # Check feature info
  if (!SeuratObject::IsMatrixEmpty(loadings)) {
    finfo <- vapply_lgl(X = list(fidx, nfeatures), FUN = is.null)
    msg <- if (all(finfo)) {
      "No feature information provided, not adding feature loadings"
    } else if (any(finfo) && !all(finfo)) {
      "Either both 'fidx' and 'nfeatures' must be supplied or both must be NULL"
    } else if (max(fidx) > nfeatures) {
      "'fidx' exceeds 'nfeatures'"
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
    varm <- if (!'varm' %in% soma_parent$names()) {
      SOMACollectionCreate(
        uri = file_path(soma_parent$uri, 'varm'),
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      )
    } else if (isTRUE(relative)) {
      SOMACollectionOpen(uri = file_path(soma_parent$uri, 'varm'), mode = 'WRITE')
    } else {
      soma_parent$varm
    }
    withCallingHandlers(
      .register_soma_object(varm, soma_parent, key = 'varm', relative = relative),
      existingKeyWarning = .maybe_muffle
    )
    on.exit(varm$close(), add = TRUE, after = FALSE)

    # Pad our feature loadings matrix
    mat <- matrix(data = NA_real_, nrow = nfeatures, ncol = ncol(loadings))
    mat[fidx, ] <- loadings

    # Write the feature loadings
    spdl::info("Adding feature loadings as {}", sQuote(ldgs))

    # Always write reductions as sparse arrays
    write_soma(
      x = mat,
      uri = ldgs,
      soma_parent = varm,
      sparse = TRUE,
      transpose = FALSE,
      key = ldgs,
      ingest_mode = ingest_mode,
      shape = dload,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  }

  return(invisible(soma_parent))
}

#' @rdname write_soma_seurat_sub
#'
#' @section Writing \code{\link[SeuratObject]{Graph}s}:
#' \pkg{Seurat} \code{\link[SeuratObject]{Graph}} objects are
#' written out as \link[tiledbsoma:SOMASparseNDArray]{sparse matrices}
#' to the \dQuote{\code{obsp}} group of a
#' \link[tiledbsoma:SOMAMeasurement]{measurement}
#'
#' @method write_soma Graph
#' @export
#'
write_soma.Graph <- function(
  x,
  uri,
  soma_parent,
  ...,
  ingest_mode = 'write',
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  relative = TRUE
) {
  check_package('SeuratObject', version = .MINIMUM_SEURAT_VERSION())
  stopifnot(
    "'soma_parent' must be a SOMAMeasurement" = inherits(
      x = soma_parent,
      what = 'SOMAMeasurement'
    ),
    "'relative' must be a single logical value" = is_scalar_logical(relative)
  )
  obsp <- if (!'obsp' %in% soma_parent$names()) {
    SOMACollectionCreate(
      uri = file_path(soma_parent$uri, 'obsp'),
      ingest_mode =  ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  } else if (isTRUE(relative)) {
    SOMACollectionOpen(uri = file_path(soma_parent$uri, 'obsp'), mode = 'WRITE')
  } else {
    soma_parent$obsp
  }
  withCallingHandlers(
    .register_soma_object(obsp, soma_parent, key = 'obsp', relative = relative),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(obsp$close(), add = TRUE, after = FALSE)

  # Find `shape` if and only if we're called from `write_soma.Seurat()`
  parents <- unique(sys.parents())
  idx <- which(vapply_lgl(
    parents,
    FUN = function(i) identical(sys.function(i), write_soma.Seurat)
  ))
  shape <- if (length(idx) == 1L) {
    get("shape", envir = sys.frame(parents[idx]))
  } else {
    NULL
  }
  if (!is.null(shape)) {
    shape <- rep_len(shape[2L], length.out = 2L)
  }

  NextMethod(
    generic = 'write_soma',
    object = x,
    uri = uri,
    soma_parent = obsp,
    sparse = TRUE,
    transpose = FALSE,
    key = basename(uri),
    ingest_mode = ingest_mode,
    shape = shape,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )

  return(invisible(soma_parent))
}

#' Write a \code{\link[SeuratObject]{Seurat}} object to a SOMA
#'
#' @inheritParams write_soma
#' @inheritParams write_soma_objects
#' @param x A \code{\link[SeuratObject]{Seurat}} object
#'
#' @inherit write_soma return
#'
#' @section Writing Cell-Level Meta Data:
#' Cell-level meta data is written out as a
#' \link[tiledbsoma:SOMADataFrame]{data frame} called \dQuote{\code{obs}} at
#' the \code{\link[tiledbsoma:SOMAExperiment]{experiment}} level
#'
#' @inherit write_soma_seurat_sub_objects sections
#'
#' @method write_soma Seurat
#' @export
#'
write_soma.Seurat <- function(
  x,
  uri,
  ...,
  ingest_mode = 'write',
  platform_config = NULL,
  tiledbsoma_ctx = NULL
) {
  check_package('SeuratObject', version = .MINIMUM_SEURAT_VERSION())
  stopifnot(
    "'uri' must be a single character value" = is.null(uri) ||
      (is_scalar_character(uri) && nzchar(uri))
  )
  ingest_mode <- match.arg(arg = ingest_mode, choices = c('write', 'resume'))
  if ('shape' %in% names(args <- rlang::dots_list(...))) {
    shape <- args$shape
    stopifnot(
      "'shape' must be a vector of two postiive integers" = is.null(shape) ||
        (rlang::is_integerish(shape, n = 2L, finite = TRUE) && all(shape > 0L))
    )
  } else {
    shape <- NULL
  }

  if (!is.null(shape) && any(shape < dim(x))) {
    stop(
      "Requested an array of shape (",
      paste(shape, collapse = ', '),
      "), but was given a Seurat object with a larger shape (",
      paste(dim(x), collapse = ', '),
      ")",
      call. = FALSE
    )
  }

  experiment <- SOMAExperimentCreate(
    uri = uri,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  on.exit(experiment$close(), add = TRUE, after = FALSE)

  # Write cell-level meta data (obs)
  spdl::info("Adding cell-level meta data")
  obs_df <- .df_index(
    x = x[[]],
    alt = 'cells',
    axis = 'obs',
    prefix = 'seurat'
  )
  obs_df[[attr(obs_df, 'index')]] <- colnames(x)
  write_soma(
    x = obs_df,
    uri = 'obs',
    soma_parent = experiment,
    key = 'obs',
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )

  # Write assays
  expms <- SOMACollectionCreate(
    uri = file_path(experiment$uri, 'ms'),
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  withCallingHandlers(
    expr = .register_soma_object(expms, soma_parent = experiment, key = 'ms'),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(expms$close(), add = TRUE, after = FALSE)

  for (measurement in SeuratObject::Assays(x)) {
    spdl::info("Adding assay {}", sQuote(measurement))
    tryCatch(
      expr = withCallingHandlers(
        .register_soma_object(
          write_soma(
            x = x[[measurement]],
            uri = measurement,
            soma_parent = expms,
            ingest_mode = ingest_mode,
            platform_config = platform_config,
            tiledbsoma_ctx = tiledbsoma_ctx
          ),
          soma_parent = expms,
          key = measurement
        ),
        existingKeyWarning = .maybe_muffle
      ),
      error = function(err) {
        if (measurement == SeuratObject::DefaultAssay(x)) {
          stop(err)
        }
        err_to_warn(err)
      }
    )
  }

  # Write dimensional reductions (obsm/varm)
  for (reduc in SeuratObject::Reductions(x)) {
    measurement <- SeuratObject::DefaultAssay(x[[reduc]])
    ms <- if (measurement %in% expms$names()) {
      SOMAMeasurementOpen(file_path(expms$uri, measurement), 'WRITE')
    } else if (SeuratObject::IsGlobal(x[[reduc]])) {
      measurement <- SeuratObject::DefaultAssay(x)
      warning(
        paste(
          strwrap(paste0(
            "Cannot find a measurement for global reduction ",
            sQuote(reduc),
            " (default assay: ",
            sQuote(SeuratObject::DefaultAssay(x[[reduc]])),
            "), adding to measurement for ",
            sQuote(measurement)
          )),
          collapse = '\n'
        ),
        call. = FALSE,
        immediate. = TRUE
      )
      SOMAMeasurementOpen(file_path(expms$uri, measurement), 'WRITE')
    } else {
      # This should never happen
      warning(
        paste(
          strwrap(paste0(
            "Cannot find a measurement for non-global reduction ",
            sQuote(reduc),
            " (default assay: ",
            sQuote(measurement),
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
      fidx <- match(
        x = rownames(loadings),
        table = rownames(x[[measurement]])
      )
      nfeatures <- nrow(x[[measurement]])
    } else {
      fidx <- nfeatures <- NULL
    }
    spdl::info("Adding dimensional reduction {}", sQuote(reduc))
    tryCatch(
      expr = write_soma(
        x = x[[reduc]],
        uri = NULL,
        soma_parent = ms,
        fidx = fidx,
        nfeatures = nfeatures,
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      ),
      error = err_to_warn
    )
    ms$close()
  }

  # Write graphs (obsp)
  for (obsp in SeuratObject::Graphs(x)) {
    measurement <- SeuratObject::DefaultAssay(x[[obsp]])
    if (!measurement %in% expms$names()) {
      warning(
        paste(
          strwrap(paste0(
            "Cannot find a measurement for graph ",
            sQuote(obsp),
            " (default assay: ",
            sQuote(measurement),
            "), skipping"
          )),
          collapse = FALSE
        ),
        call. = FALSE,
        immediate. = TRUE
      )
      next
    }
    ms <- SOMAMeasurementOpen(file_path(expms$uri, measurement))
    spdl::info("Adding graph {}", sQuote(obsp))
    tryCatch(
      expr = write_soma(
        x = x[[obsp]],
        uri = obsp,
        soma_parent = ms,
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      ),
      error = err_to_warn
    )
    ms$close()
  }

  # TODO: Write images
  if (length(SeuratObject::Images(x))) {
    warning(
      "Spatially resolved data cannot be written to SOMAs at this time",
      call. = FALSE,
      immediate. = TRUE
    )
  }

  # Add extra Seurat data
  expuns <- SOMACollectionCreate(
    uri = file_path(experiment$uri, 'uns'),
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  withCallingHandlers(
    expr = .register_soma_object(expuns, soma_parent = experiment, key = 'uns'),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(expuns$close(), add = TRUE, after = FALSE)

  # Write command logs
  for (cmd in SeuratObject::Command(x)) {
    spdl::info("Adding command log '{}'", cmd)
    write_soma(
      x = x[[cmd]],
      uri = cmd,
      soma_parent = expuns,
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  }

  return(experiment$uri)
}

#' @rdname write_soma_seurat_sub
#'
#' @section Writing \code{\link[SeuratObject]{SeuratCommand}s}:
#' \pkg{Seurat} \link[SeuratObject:SeuratCommand]{command logs} are written out
#' as \link[tiledbsoma:SOMADataFrame]{data frames} to the
#' \dQuote{\code{seurat_commands}} group of a
#' \link[tiledbsoma:SOMACollection]{collection}
#'
#' @method write_soma SeuratCommand
#' @export
#'
write_soma.SeuratCommand <- function(
  x,
  uri = NULL,
  soma_parent,
  ...,
  ingest_mode = 'write',
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  relative = TRUE
) {
  check_package('SeuratObject', version = .MINIMUM_SEURAT_VERSION())
  check_package('jsonlite')
  stopifnot(
    "'uri' must be a single character value" = is.null(uri) ||
      (is_scalar_character(uri) && nzchar(uri)),
    "'soma_parent' must be a SOMACollection" = inherits(soma_parent, what = 'SOMACollection'),
    "'relative' must be a single logical value" = is_scalar_logical(relative)
  )

  key <- 'seurat_commands'
  uri <- uri %||% methods::slot(x, name = 'name')

  # Create a group for command logs
  logs_uri <- .check_soma_uri(key, soma_parent = soma_parent, relative = relative)
  logs <- if (!key %in% soma_parent$names()) {
    spdl::info("Creating a group for command logs")
    logs <- SOMACollectionCreate(
      uri = logs_uri,
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
    soma_parent$add_new_collection(logs, key)
    logs
  } else {
    logs <- soma_parent$get(key)
    if (!inherits(logs, 'SOMACollection')) {
      stop(
        "Existing ",
        class(logs)[1L],
        " named ",
        sQuote(key),
        " found, expected a SOMACollection",
        call. = FALSE
      )
    }
    if (isTRUE(relative)) {
      logs$close()
      logs <- SOMACollectionOpen(logs_uri, mode = 'WRITE')
    }
    spdl::info("Found existing group for command logs")
    logs$reopen("WRITE")
    logs
  }
  on.exit(logs$close(), add = TRUE, after = FALSE)

  # Encode parameters
  spdl::info("Encoding parameters in the command log")
  xlist <- as.list(x, complete = TRUE)
  for (i in names(xlist)) {
    # Timestamp -> JSON defaults to:
    # - timestamp -> string
    # - string -> JSON
    # which is lossy. Instead, do
    # - timestamp -> numeric
    # - numeric -> hex double precision (`sprintf("%a")`)
    # - hex double precision -> JSON
    # for lossless timestamp encoding in JSON
    if (i == 'time.stamp') {
      ts <- sapply(
        unclass(as.POSIXlt(
          xlist[[i]],
          tz = attr(xlist[[i]], 'tzone', exact = TRUE) %||% Sys.timezone()
        )),
        .encode_as_char,
        simplify = FALSE,
        USE.NAMES = TRUE
      )
      xlist[[i]] <- as.character(jsonlite::toJSON(ts, auto_unbox = TRUE))
    }
    # Encode numerics/doubles as hex double precision for lossless encoding
    xlist[[i]] <- .encode_as_char(xlist[[i]])
  }

  # Encode as JSON
  spdl::info("Encoding command log as JSON")
  enc <- as.character(jsonlite::toJSON(
    xlist,
    null = 'null',
    auto_unbox = TRUE
  ))

  # Write out and return
  sdf <- write_soma(
    x = enc,
    uri = uri,
    soma_parent = logs,
    key = basename(uri),
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    relative = relative
  )

  return(invisible(soma_parent))
}
