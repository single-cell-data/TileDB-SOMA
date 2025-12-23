#' Convert a \pkg{Seurat} Sub-Object to a SOMA Object, returned opened for write
#'
#' Various helpers to write \pkg{Seurat} sub-objects to SOMA objects..
#'
#' @inheritParams write_soma_objects
#' @param x A \pkg{Seurat} sub-object (eg. an
#' \code{\link[SeuratObject]{Assay}}, a \code{\link[SeuratObject]{DimReduc}},
#' or a \code{\link[SeuratObject]{Graph}}).
#' @param soma_parent The parent \link[tiledbsoma:SOMACollection]{collection};
#' for the \code{DimReduc} and \code{Graph} methods, this \strong{must} be a
#' \link[tiledbsoma:SOMAMeasurement]{measurement} for the assay \code{x}
#' was generated from.
#'
#' @name write_soma_seurat_sub_objects
#' @rdname write_soma_seurat_sub
#'
#' @keywords internal
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE) && requireNamespace("SeuratObject", quietly = TRUE)
#' uri <- withr::local_tempfile(pattern = "seurat-sub")
#'
#' data("pbmc_small", package = "SeuratObject")
#' suppressWarnings(pbmc_small <- SeuratObject::UpdateSeuratObject(pbmc_small))
#'
#' col <- SOMACollectionCreate(uri)
#'
#' # Write a v3 Assay
#' (assay <- pbmc_small[["RNA"]])
#' (ms <- write_soma(assay, "RNA", soma_parent = col))
#'
#' # Write a v5 Assay
#' (assay5 <- methods::as(pbmc_small[["RNA"]], "Assay5"))
#' (ms5 <- write_soma(assay5, "RNA5", soma_parent = col))
#'
#' ms5$close()
#'
#' # Write a dimensional reduction
#' (tsne <- pbmc_small[["tsne"]])
#' write_soma(tsne, soma_parent = ms)
#' ms$obsm
#'
#' # Write a Seurat Graph
#' (snn <- pbmc_small[["RNA_snn"]])
#' write_soma(snn, "snn", soma_parent = ms)
#' ms$obsp
#'
#' # Write a Seurat command log
#' (cmd <- pbmc_small[["NormalizeData.RNA"]])
#' write_soma(cmd, "NormalizeData.RNA", soma_parent = col)
#' (logs <- col$get("seurat_commands"))
#' logs$get("NormalizeData.RNA")
#'
#' col$close()
#'
NULL

#' @return \code{Assay} and \code{Assay5} methods: a
#' \code{\link{SOMAMeasurement}} with the data from \code{x},
#' returned opened for write.
#'
#' @rdname write_soma_seurat_sub
#'
#' @section Writing v3 \code{\link[SeuratObject]{Assay}s}:
#' \pkg{Seurat} \code{\link[SeuratObject]{Assay}} objects are written out as
#' individual \link[tiledbsoma:SOMAMeasurement]{measurements}:
#' \itemize{
#'  \item the \dQuote{\code{data}} matrix is written out as a
#'   \link[tiledbsoma:SOMASparseNDArray]{sparse array} called
#'   \dQuote{\code{data}} within the \dQuote{\code{X}} group.
#'  \item the \dQuote{\code{counts}} matrix, if not
#'   \link[SeuratObject:IsMatrixEmpty]{empty}, is written out as a
#'   \link[tiledbsoma:SOMASparseNDArray]{sparse array} called
#'   \dQuote{\code{counts}} within the \dQuote{\code{X}} group.
#'  \item the \dQuote{\code{scale.data}} matrix, if not
#'   \link[SeuratObject:IsMatrixEmpty]{empty}, is written out as a
#'   \link[tiledbsoma:SOMASparseNDArray]{sparse array} called
#'   \dQuote{\code{scale_data}} within the \dQuote{\code{X}} group.
#'  \item feature-level metadata is written out as a
#'   \link[tiledbsoma:SOMADataFrame]{data frame} called \dQuote{\code{var}}.
#' }
#' Expression matrices are transposed (cells as rows) prior to writing. All
#' other slots, including results from extended assays (eg. \code{SCTAssay},
#' \code{ChromatinAssay}) are lost.
#'
#' \subsection{Performance Considerations}{
#' Ingestion of very large dense layers, such as \code{scale.data}, can be
#' memory intensive. For better performance, users can remove these layers
#' prior to ingestion and regenerate them after export, or ingest them
#' separately as dense arrays for those who need to persist the exact matrix
#'
#' \preformatted{
#' # Using SeuratObject v5 syntax on a v3 `Assay`
#' # Cache the layer for separate ingestion, skip if planning to regenerate
#' mat <- object[["ASSAY"]]$scale.data
#'
#' # Remove the `scale.data` layer
#' object[["ASSAY"]]$scale.data <- NULL
#'
#' # Ingest the smaller object
#' uri <- write_soma(object, "/path/to/soma")
#'
#' # Ingest the `scale.data` layer densely; needed only if persistence
#' # of the data is paramount
#' # Pad the `scale.data` layer so that its soma join IDs match the experiment
#' padded <- matrix(
#'   data = vector("numeric", length = prod(dim(object[["ASSAY"]]))),
#'   nrow = nrow(object[["ASSAY"]]),
#'   ncol = ncol(object[["ASSAY"]])
#' )
#' rowidx <- match(rownames(mat), rownames(object[["ASSAY"]]))
#' colidx <- match(colnames(mat), colnames(object[["ASSAY"]]))
#' padded[rowidx, colidx] <- mat
#'
#' # Use `write_soma()` to ingest densely and register it within the `uns`
#' # collection; this may need to be created manually if the original
#' # object does not contain command logs
#' exp <- SOMAExperimentOpen(uri, "WRITE")
#' if (!match("uns", exp$names(), nomatch = 0L)) {
#'   # For `tiledb://` URIs, set the URI for the new collection manually rather
#'   # than relying on `file.path()`
#'   uns <- SOMACollectionCreate(file.path(exp$uri, "uns"))
#'   exp$add_new_collection(uns, "uns")
#' }
#' arr <- write_soma(
#'   padded,
#'   "scale_data",
#'   soma_parent = exp$get("uns"),
#'   sparse = FALSE,
#'   key = "scale_data"
#' )
#' arr$close()
#' exp$close()
#' }
#' Please note that dense arrays cannot be read in using the
#' \code{\link{SOMAExperimentAxisQuery}} mechanism; use
#' \code{\link{SOMADenseNDArray}$read_dense_matrix}, remembering to transpose
#' before adding back to a \code{Seurat} object
#' }
#'
#' @method write_soma Assay
#' @export
#'
write_soma.Assay <- .write_seurat_assay

#' @rdname write_soma_seurat_sub
#'
#' @section Writing v5 \code{Assays}:
#' \pkg{Seurat} v5 \code{\link[SeuratObject:Assay5]{Assays}s} are written
#' out as individual \link[tiledbsoma:SOMAMeasurement]{measurements}:
#' \itemize{
#'  \item the layer matrices are written out as
#'   \link[tiledbsoma:SOMASparseNDArray]{sparse arrays} within the
#'   \dQuote{\code{X}} group.
#'  \item feature-level metadata is written out as a
#'   \link[tiledbsoma:SOMADataFrame]{data frame} called \dQuote{\code{var}}.
#' }
#' Expression matrices are transposed (cells as rows) prior to writing. All
#' other slots, including results from extended assays (eg. \code{SCTAssay},
#' \code{ChromatinAssay}) are lost.\cr
#' The following bits of metadata are written in various parts of the measurement
#' \itemize{
#'  \item \dQuote{\code{soma_ecosystem_seurat_assay_version}}: written at the
#'   measurement level; indicates the Seurat assay version.
#'   Set to \dQuote{\code{v5}}.
#'  \item \dQuote{\code{soma_ecosystem_seurat_v5_default_layers}}: written at
#'   the \dQuote{\code{X}} group level; indicates the
#'   \link[SeuratObject:DefaultLayer]{default layers}.
#'  \item \dQuote{\code{soma_ecosystem_seurat_v5_ragged}}: written at the
#'   \dQuote{\code{X/<layer>}} array level; with a value of
#'   \dQuote{\code{ragged}}, indicates whether or not the layer is ragged.
#'  \item \dQuote{\code{soma_r_type_hint}}: written at the
#'   \dQuote{\code{X/<layer>}} array level; indicates the \R class and
#'   defining package (for S4 classes) of the original layer.
#' }
#'
#' @method write_soma Assay5
#' @export
#'
write_soma.Assay5 <- .write_seurat_assay

#' @param fidx An integer vector describing the location of features in
#' \code{SeuratObject::Loadings(x)} with relation to \code{soma_parent}
#' (eg. \code{match(rownames(Loadings(x)), rownames(assay))}).
#' @param nfeatures The number of features present in \code{soma_parent}
#' (eg. \code{nrow(assay)}).
#'
#' @return \code{DimReduc} and \code{Graph} methods: invisibly returns
#' \code{soma_parent}, opened for write, with the values of \code{x}
#' added to it.
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
#'   \dQuote{\code{obsm}} group.
#'  \item feature loadings, if not \link[SeuratObject:IsMatrixEmpty]{empty},
#'   are written out as a \link[tiledbsoma:SOMASparseNDArray]{sparse matrix} in
#'   the \dQuote{\code{varm}} groups; loadings are padded with \code{NAs}
#'   to include all features.
#' }
#' Dimensional reduction names are translated to AnnData-style names (eg.
#' \dQuote{\code{pca}} becomes \code{X_pca} for embeddings and
#' \dQuote{\code{PCs}} for loadings). All other slots, including projected
#' feature loadings and jackstraw information, are lost.
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
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL,
  relative = TRUE
) {
  check_package("SeuratObject", version = .MINIMUM_SEURAT_VERSION())
  stopifnot(
    "'uri' must be NULL" = is.null(uri),
    "'soma_parent' must be a SOMAMeasurement" = inherits(
      x = soma_parent,
      what = "SOMAMeasurement"
    ),
    "'fidx' must be a positive integer vector" = is.null(fidx) ||
      (rlang::is_integerish(fidx, finite = TRUE) && all(fidx > 0L)),
    "'nfeatures' must be a single positive integer" = is.null(nfeatures) ||
      (rlang::is_integerish(nfeatures, n = 1L, finite = TRUE) &&
        nfeatures > 0L),
    "'relative' must be a single logical value" = is_scalar_logical(relative)
  )

  key <- tolower(gsub(
    pattern = "_$",
    replacement = "",
    x = SeuratObject::Key(x)
  ))
  key <- switch(EXPR = key, pc = "pca", ic = "ica", key)

  # Find `shape` if and only if we're called from `write_soma.Seurat()`
  parents <- unique(sys.parents())
  idx <- which(vapply_lgl(
    X = parents,
    FUN = function(i) {
      return(identical(sys.function(i), write_soma.Seurat))
    }
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
  obsm <- if (!"obsm" %in% soma_parent$names()) {
    SOMACollectionCreate(
      uri = file_path(soma_parent$uri, "obsm"),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      context = context
    )
  } else if (isTRUE(relative)) {
    SOMACollectionOpen(uri = file_path(soma_parent$uri, "obsm"), mode = "WRITE")
  } else {
    soma_parent$obsm
  }
  withCallingHandlers(
    .register_soma_object(obsm, soma_parent, key = "obsm", relative = relative),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(obsm$close(), add = TRUE, after = FALSE)

  embed <- paste0("X_", key)
  soma_info(sprintf("Adding embeddings as %s", sQuote(embed)))

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
    tiledbsoma_ctx = tiledbsoma_ctx,
    context = context
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
      ""
    }
    if (nzchar(msg)) {
      warning(
        paste(
          strwrap(paste0(msg, ", not adding feature loadings")),
          collapse = "\n"
        ),
        call. = FALSE,
        immediate. = TRUE
      )
      loadings <- methods::new("matrix")
    }
  }

  # Write feature loadings
  if (!SeuratObject::IsMatrixEmpty(loadings)) {
    ldgs <- switch(
      EXPR = key,
      pca = "PCs",
      ica = "ICs",
      paste0(toupper(key), "s")
    )

    # Create a group for `varm`
    varm <- if (!"varm" %in% soma_parent$names()) {
      SOMACollectionCreate(
        uri = file_path(soma_parent$uri, "varm"),
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx,
        context = context
      )
    } else if (isTRUE(relative)) {
      SOMACollectionOpen(
        uri = file_path(soma_parent$uri, "varm"),
        mode = "WRITE"
      )
    } else {
      soma_parent$varm
    }
    withCallingHandlers(
      .register_soma_object(
        varm,
        soma_parent,
        key = "varm",
        relative = relative
      ),
      existingKeyWarning = .maybe_muffle
    )
    on.exit(varm$close(), add = TRUE, after = FALSE)

    # Pad our feature loadings matrix
    mat <- matrix(data = NA_real_, nrow = nfeatures, ncol = ncol(loadings))
    mat[fidx, ] <- loadings

    # Write the feature loadings
    soma_info(sprintf("Adding feature loadings as %s", sQuote(ldgs)))

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
      tiledbsoma_ctx = tiledbsoma_ctx,
      context = context
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
#' \link[tiledbsoma:SOMAMeasurement]{measurement}.
#'
#' @method write_soma Graph
#' @export
#'
write_soma.Graph <- function(
  x,
  uri,
  soma_parent,
  ...,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL,
  relative = TRUE
) {
  check_package("SeuratObject", version = .MINIMUM_SEURAT_VERSION())
  stopifnot(
    "'soma_parent' must be a SOMAMeasurement" = inherits(
      x = soma_parent,
      what = "SOMAMeasurement"
    ),
    "'relative' must be a single logical value" = is_scalar_logical(relative)
  )
  obsp <- if (!"obsp" %in% soma_parent$names()) {
    SOMACollectionCreate(
      uri = file_path(soma_parent$uri, "obsp"),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      context = context
    )
  } else if (isTRUE(relative)) {
    SOMACollectionOpen(uri = file_path(soma_parent$uri, "obsp"), mode = "WRITE")
  } else {
    soma_parent$obsp
  }
  withCallingHandlers(
    .register_soma_object(obsp, soma_parent, key = "obsp", relative = relative),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(obsp$close(), add = TRUE, after = FALSE)

  # Find `shape` if and only if we're called from `write_soma.Seurat()`
  parents <- unique(sys.parents())
  idx <- which(vapply_lgl(
    X = parents,
    FUN = function(i) {
      return(identical(sys.function(i), write_soma.Seurat))
    }
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
    generic = "write_soma",
    object = x,
    uri = uri,
    soma_parent = obsp,
    sparse = TRUE,
    type = arrow::infer_type(methods::slot(x, "x")),
    transpose = FALSE,
    key = basename(uri),
    ingest_mode = ingest_mode,
    shape = shape,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    context = context
  )

  return(invisible(soma_parent))
}

#' Write a \code{\link[SeuratObject]{Seurat}} object to a SOMA
#'
#' @inheritParams write_soma
#' @inheritParams write_soma_objects
#' @param x A \code{\link[SeuratObject]{Seurat}} object.
#'
#' @inherit write_soma return
#'
#' @section Writing Cell-Level Metadata:
#' Cell-level metadata is written out as a
#' \link[tiledbsoma:SOMADataFrame]{data frame} called \dQuote{\code{obs}} at
#' the \code{\link[tiledbsoma:SOMAExperiment]{experiment}} level.
#'
#' @inherit write_soma_seurat_sub_objects sections
#'
#' @method write_soma Seurat
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE) && requireNamespace("SeuratObject", quietly = TRUE)
#' \donttest{
#' uri <- withr::local_tempfile(pattern = "pbmc-small")
#'
#' data("pbmc_small", package = "SeuratObject")
#' suppressWarnings(pbmc_small <- SeuratObject::UpdateSeuratObject(pbmc_small))
#'
#' uri <- write_soma(pbmc_small, uri)
#'
#' (exp <- SOMAExperimentOpen(uri))
#' exp$obs
#' exp$get("uns")$get("seurat_commands")$names()
#' (ms <- exp$ms$get("RNA"))
#' ms$var
#' ms$X$names()
#' ms$obsm$names()
#' ms$varm$names()
#' ms$obsp$names()
#'
#' exp$close()
#' }
#'
write_soma.Seurat <- function(
  x,
  uri,
  ...,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL
) {
  context = get_soma_context(context, tiledbsoma_ctx, what="Seurat.write_soma(tiledbsoma_ctx)")
  # Allow writing `soma_` prefixed columns to SOMADataFrames
  # (normally disallowed as a reserved prefix)
  op <- options(tiledbsoma.write_soma.internal = TRUE)
  on.exit(options(op), add = TRUE, after = FALSE)
  check_package("SeuratObject", version = .MINIMUM_SEURAT_VERSION())
  stopifnot(
    "'uri' must be a single character value" = is.null(uri) ||
      (is_scalar_character(uri) && nzchar(uri))
  )
  ingest_mode <- match.arg(arg = ingest_mode, choices = c("write", "resume"))
  if ("shape" %in% names(args <- rlang::dots_list(...))) {
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
      paste(shape, collapse = ", "),
      "), but was given a Seurat object with a larger shape (",
      paste(dim(x), collapse = ", "),
      ")",
      call. = FALSE
    )
  }

  experiment <- SOMAExperimentCreate(
    uri = uri,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    context = context
  )
  on.exit(experiment$close(), add = TRUE, after = FALSE)

  # Prepare cell-level metadata (obs)
  obs_df <- .df_index(x = x[[]], alt = "cells", axis = "obs", prefix = "seurat")
  obs_df[[attr(obs_df, 'index')]] <- colnames(x)

  # Write assays
  expms <- SOMACollectionCreate(
    uri = file_path(experiment$uri, "ms"),
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    context = context
  )
  withCallingHandlers(
    expr = .register_soma_object(expms, soma_parent = experiment, key = "ms"),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(expms$close(), add = TRUE, after = FALSE)

  for (measurement in SeuratObject::Assays(x)) {
    soma_info(sprintf("Adding assay %s", sQuote(measurement)))
    tryCatch(
      expr = withCallingHandlers(
        .register_soma_object(
          write_soma(
            x = x[[measurement]],
            uri = measurement,
            soma_parent = expms,
            ingest_mode = ingest_mode,
            platform_config = platform_config,
            tiledbsoma_ctx = tiledbsoma_ctx,
            context = context
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
    if (inherits(x[[measurement]], what = 'StdAssay')) {
      key <- .assay_obs_hint(measurement)
      obs_df[[key]] <- colnames(x) %in% colnames(x[[measurement]])
    }
  }

  # Write cell-level metadata (obs)
  soma_info("Adding cell-level metadata")
  write_soma(
    x = obs_df,
    uri = 'obs',
    soma_parent = experiment,
    key = 'obs',
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    context = context
  )

  # Write dimensional reductions (obsm/varm)
  expms$reopen("WRITE")
  for (reduc in SeuratObject::Reductions(x)) {
    measurement <- SeuratObject::DefaultAssay(x[[reduc]])
    ms <- if (measurement %in% expms$names()) {
      SOMAMeasurementOpen(file_path(expms$uri, measurement), "WRITE")
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
          collapse = "\n"
        ),
        call. = FALSE,
        immediate. = TRUE
      )
      SOMAMeasurementOpen(file_path(expms$uri, measurement), "WRITE")
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
          collapse = "\n"
        ),
        call. = FALSE,
        immediate. = TRUE
      )
      next
    }
    loadings <- SeuratObject::Loadings(x[[reduc]])
    if (!SeuratObject::IsMatrixEmpty(loadings)) {
      fidx <- match(x = rownames(loadings), table = rownames(x[[measurement]]))
      nfeatures <- nrow(x[[measurement]])
    } else {
      fidx <- nfeatures <- NULL
    }
    soma_info(sprintf("Adding dimensional reduction %s", sQuote(reduc)))
    tryCatch(
      expr = write_soma(
        x = x[[reduc]],
        uri = NULL,
        soma_parent = ms,
        fidx = fidx,
        nfeatures = nfeatures,
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx,
        context = context
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
    soma_info(sprintf("Adding graph %s", sQuote(obsp)))
    tryCatch(
      expr = write_soma(
        x = x[[obsp]],
        uri = obsp,
        soma_parent = ms,
        ingest_mode = ingest_mode,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx,
        context = context
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
    uri = file_path(experiment$uri, "uns"),
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    context = context
  )
  withCallingHandlers(
    expr = .register_soma_object(expuns, soma_parent = experiment, key = "uns"),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(expuns$close(), add = TRUE, after = FALSE)

  # Write command logs
  for (cmd in SeuratObject::Command(x)) {
    soma_info(sprintf("Adding command log '%s'", cmd))
    write_soma(
      x = x[[cmd]],
      uri = cmd,
      soma_parent = expuns,
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      context = context
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
#' \link[tiledbsoma:SOMACollection]{collection}.
#'
#' @method write_soma SeuratCommand
#' @export
#'
write_soma.SeuratCommand <- function(
  x,
  uri = NULL,
  soma_parent,
  ...,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL,
  relative = TRUE
) {
  check_package("SeuratObject", version = .MINIMUM_SEURAT_VERSION())
  check_package("jsonlite")
  stopifnot(
    "'uri' must be a single character value" = is.null(uri) ||
      (is_scalar_character(uri) && nzchar(uri)),
    "'soma_parent' must be a SOMACollection" = inherits(
      soma_parent,
      what = "SOMACollection"
    ),
    "'relative' must be a single logical value" = is_scalar_logical(relative)
  )

  key <- "seurat_commands"
  uri <- uri %||% methods::slot(x, name = "name")

  # Create a group for command logs
  logs_uri <- .check_soma_uri(
    key,
    soma_parent = soma_parent,
    relative = relative
  )
  logs <- if (!key %in% soma_parent$names()) {
    soma_info("Creating a group for command logs")
    logs <- SOMACollectionCreate(
      uri = logs_uri,
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      context = context
    )
    soma_parent$add_new_collection(logs, key)
    logs
  } else {
    logs <- soma_parent$get(key)
    if (!inherits(logs, "SOMACollection")) {
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
      logs <- SOMACollectionOpen(logs_uri, mode = "WRITE")
    }
    soma_info("Found existing group for command logs")
    logs$reopen("WRITE")
    logs
  }
  on.exit(logs$close(), add = TRUE, after = FALSE)

  # Encode parameters
  soma_info("Encoding parameters in the command log")
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
    if (i == "time.stamp") {
      ts <- sapply(
        unclass(as.POSIXlt(
          xlist[[i]],
          tz = attr(xlist[[i]], "tzone", exact = TRUE) %||% Sys.timezone()
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
  soma_info("Encoding command log as JSON")
  enc <- as.character(jsonlite::toJSON(xlist, null = "null", auto_unbox = TRUE))

  # Write out and return
  sdf <- write_soma(
    x = enc,
    uri = uri,
    soma_parent = logs,
    key = basename(uri),
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    context = context,
    relative = relative
  )

  return(invisible(soma_parent))
}
