#' @name write_soma_objects
#' @rdname write_soma_objects
#'
#' @method write_soma DataFrame
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE) && requireNamespace("SeuratObject", quietly = TRUE) && requireNamespace("S4Vectors", quietly = TRUE)
#' # Write a Bioconductor S4 DataFrame object to a SOMA
#' uri <- withr::local_tempfile(pattern = "s4-data-frame")
#' data("pbmc_small", package = "SeuratObject")
#' obs <- suppressWarnings(SeuratObject::UpdateSeuratObject(pbmc_small))[[]]
#' head(obs <- as(obs, "DataFrame"))
#'
#' (sdf <- write_soma(obs, uri, soma_parent = NULL, relative = FALSE))
#'
#' sdf$close()
#'
write_soma.DataFrame <- function(
  x,
  uri,
  soma_parent,
  df_index = NULL,
  index_column_names = "soma_joinid",
  ...,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  soma_context = NULL,
  relative = TRUE
) {
  # Check for compound non-atomic/factor types
  for (i in names(x)) {
    if (!(is.atomic(x[[i]]) || is.factor(x[[i]]))) {
      stop("All columns in DataFrames must be atomic or factors", call. = FALSE)
    }
  }
  index <- attr(x, which = "index")
  x <- suppressWarnings(as.data.frame(x), classes = "deprecatedWarning")
  attr(x, which = "index") <- index
  return(write_soma(
    x = x,
    uri = uri,
    soma_parent = soma_parent,
    df_index = df_index,
    index_column_names = index_column_names,
    ...,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    soma_context = soma_context,
    relative = relative
  ))
}

#' @rdname write_soma_objects
#'
#' @method write_soma Hits
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE) && requireNamespace("S4Vectors", quietly = TRUE)
#' # Write a Bioconductor SelfHits object to a SOMA
#' uri <- withr::local_tempfile(pattern = "hits")
#' (hits <- S4Vectors::SelfHits(
#'   c(2, 3, 3, 3, 3, 3, 4, 4, 4),
#'   c(4, 3, 2:4, 2, 2:3, 2),
#'   4,
#'   x = stats::rnorm(9L)
#' ))
#'
#' (arr <- write_soma(hits, uri, soma_parent = NULL, relative = FALSE))
#'
#' arr$close()
#'
write_soma.Hits <- function(
  x,
  uri,
  soma_parent,
  sparse = TRUE,
  type = NULL,
  transpose = FALSE,
  ...,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  soma_context = NULL,
  relative = TRUE
) {
  return(write_soma(
    x = .hits_to_mat(x),
    uri = uri,
    soma_parent = soma_parent,
    sparse = sparse,
    type = type,
    transpose = transpose,
    ...,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    soma_context = soma_context,
    relative = relative
  ))
}

#' Write a \code{\link[SingleCellExperiment:SingleCellExperiment-class]{SingleCellExperiment}}
#' object to a SOMA
#'
#' @inheritParams write_soma
#' @inheritParams write_soma_objects
#' @param ms_name Name for resulting measurement; defaults to
#' \code{\link[SingleCellExperiment]{mainExpName}(x)}.
#'
#' @inherit write_soma.SummarizedExperiment return sections
#'
#' @section Writing Reduced Dimensions:
#' Reduced dimensions are written out as
#' \link[tiledbsoma:SOMASparseNDArray]{sparse matrices} within the \code{obsm}
#' group of \code{\link[tiledbsoma:SOMAMeasurement]{measurement}}
#' named \code{ms_name}.
#'
#' @section Writing Column Pairs:
#' Column-wise relationship matrices are written out as
#' \link[tiledbsoma:SOMASparseNDArray]{sparse matrices} within the
#' \code{obsp} group of \code{\link[tiledbsoma:SOMAMeasurement]{measurement}}
#' named \code{ms_name}.
#'
#' @section Writing Row Pairs:
#' Row-wise relationship matrices are written out as
#' \link[tiledbsoma:SOMASparseNDArray]{sparse matrices} within the
#' \code{varp} group of \code{\link[tiledbsoma:SOMAMeasurement]{measurement}}
#' named \code{ms_name}.
#'
#' @method write_soma SingleCellExperiment
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE) && requireNamespace("SingleCellExperiment", quietly = TRUE)
#' \donttest{
#' uri <- withr::local_tempfile(pattern = "single-cell-experiment")
#'
#' mat <- abs(Matrix::rsparsematrix(
#'   230L,
#'   80L,
#'   0.3,
#'   dimnames = list(paste0("feature_", seq_len(230)), paste0("cell_", seq_len(80)))
#' ))
#' (sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(counts = mat, logcounts = log2(mat + 1L)),
#'   reducedDims = list(
#'     pca = matrix(stats::runif(80 * 5L), nrow = 80),
#'     tsne = matrix(stats::rnorm(80 * 2L), nrow = 80)
#'   ),
#'   mainExpName = "RNA"
#' ))
#'
#' uri <- write_soma(sce, uri)
#'
#' (exp <- SOMAExperimentOpen(uri))
#' exp$close()
#' }
#'
write_soma.SingleCellExperiment <- function(
  x,
  uri,
  ms_name = NULL,
  ...,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  soma_context = NULL
) {
  check_package("SingleCellExperiment", version = .MINIMUM_SCE_VERSION())
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
  ms_name <- ms_name %||% SingleCellExperiment::mainExpName(x)

  uri <- NextMethod(
    "write_soma",
    x,
    uri = uri,
    ms_name = ms_name,
    ...,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    soma_context = soma_context
  )
  experiment <- SOMAExperimentOpen(
    uri = uri,
    mode = "WRITE",
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    soma_context = soma_context
  )
  on.exit(expr = experiment$close(), add = TRUE, after = FALSE)

  ms <- SOMAMeasurementOpen(
    file_path(experiment$uri, "ms", ms_name),
    mode = "WRITE"
  )
  on.exit(ms$close(), add = TRUE, after = FALSE)

  # Write reduced dimensions
  soma_info("Adding reduced dimensions")
  obsm <- if (!"obsm" %in% ms$names()) {
    SOMACollectionCreate(
      uri = file_path(ms$uri, "obsm"),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      soma_context = soma_context
    )
  } else {
    SOMACollectionOpen(file_path(ms$uri, "obsm"), mode = "WRITE")
  }
  withCallingHandlers(
    .register_soma_object(obsm, soma_parent = ms, key = "obsm"),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(obsm$close(), add = TRUE, after = FALSE)

  for (rd in SingleCellExperiment::reducedDimNames(x)) {
    soma_info(sprintf("Adding reduced dimension %s", rd))
    write_soma(
      x = SingleCellExperiment::reducedDim(x, rd),
      uri = rd,
      soma_parent = obsm,
      sparse = TRUE,
      key = rd,
      ingest_mode = ingest_mode,
      shape = if (is.null(shape)) {
        NULL
      } else {
        c(shape[2L], ncol(SingleCellExperiment::reducedDim(x, rd)))
      },
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      soma_context = soma_context
    )
  }

  # Write nearest-neighbor graphs
  obsp <- if (!"obsp" %in% ms$names()) {
    SOMACollectionCreate(
      uri = file_path(ms$uri, "obsp"),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      soma_context = soma_context
    )
  } else {
    SOMACollectionOpen(file_path(ms$uri, "obsp"), mode = "WRITE")
  }
  withCallingHandlers(
    .register_soma_object(obsp, soma_parent = ms, key = "obsp"),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(obsp$close(), add = TRUE, after = FALSE)

  for (cp in SingleCellExperiment::colPairNames(x)) {
    soma_info(sprintf("Adding colPair %s", cp))
    write_soma(
      x = SingleCellExperiment::colPair(x, cp),
      uri = cp,
      soma_parent = obsp,
      sparse = TRUE,
      key = cp,
      ingest_mode = ingest_mode,
      shape = if (is.null(shape)) {
        NULL
      } else {
        rep_len(shape[2L], length.out = 2L)
      },
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      soma_context = soma_context
    )
  }

  # Write coexpression networks
  varp <- if (!"varp" %in% ms$names()) {
    SOMACollectionCreate(
      uri = file_path(ms$uri, "varp"),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      soma_context = soma_context
    )
  } else {
    SOMACollectionOpen(file_path(ms$uri, "varp"), mode = "WRITE")
  }
  withCallingHandlers(
    .register_soma_object(varp, soma_parent = ms, key = "varp"),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(varp$close(), add = TRUE, after = FALSE)

  for (rp in SingleCellExperiment::rowPairNames(x)) {
    soma_info(sprintf("Adding rowPair %s", rp))
    write_soma(
      x = SingleCellExperiment::rowPair(x, rp),
      uri = rp,
      soma_parent = varp,
      sparse = TRUE,
      key = rp,
      ingest_mode = ingest_mode,
      shape = if (is.null(shape)) {
        NULL
      } else {
        rep_len(shape[1L], length.out = 2L)
      },
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      soma_context = soma_context
    )
  }

  # TODO: Add alternate experiments
  return(experiment$uri)
}

#' Write a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object to a SOMA
#'
#' @inheritParams write_soma
#' @inheritParams write_soma_objects
#' @param ms_name Name for resulting measurement.
#'
#' @inherit write_soma return
#'
#' @section Writing \code{colData}:
#' \code{colData} is written out as a
#' \link[tiledbsoma:SOMADataFrame]{data frame} called \dQuote{\code{obs}} at
#' the \code{\link[tiledbsoma:SOMAExperiment]{experiment}} level.
#'
#' @section Writing Assay Matrices:
#' Each \link[SummarizedExperiment:assay]{assay matrix} is written out as a
#' \link[tiledbsoma:SOMASparseNDArray]{sparse matrix} within the \code{X} group
#' of \code{\link[tiledbsoma:SOMAMeasurement]{measurement}} named
#' \code{ms_name}. Names for assay matrices within \code{X} are taken from the
#' \link[SummarizedExperiment:assayNames]{assay names}. Assay matrices are
#' transposed (samples as rows) prior to writing.
#'
#' @section Writing \code{rowData}:
#' \code{rowData} is written out as a
#' \link[tiledbsoma:SOMADataFrame]{data frame} called \dQuote{\code{var}} at
#' the \code{\link[tiledbsoma:SOMAMeasurement]{measurement}} level.
#'
#' @method write_soma SummarizedExperiment
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE) && requireNamespace("SummarizedExperiment", quietly = TRUE)
#' \donttest{
#' uri <- withr::local_tempfile(pattern = "summarized-experiment")
#'
#' mat <- abs(Matrix::rsparsematrix(
#'   230L,
#'   80L,
#'   0.3,
#'   dimnames = list(paste0("feature_", seq_len(230)), paste0("cell_", seq_len(80)))
#' ))
#' (se <- SummarizedExperiment::SummarizedExperiment(list(counts = mat, logcounts = log2(mat + 1L))))
#'
#' uri <- write_soma(se, uri, ms_name = "RNA")
#'
#' (exp <- SOMAExperimentOpen(uri))
#' exp$obs
#' (ms <- exp$ms$get("RNA"))
#' ms$var
#' ms$X$names()
#'
#' exp$close()
#' }
#'
write_soma.SummarizedExperiment <- function(
  x,
  uri,
  ms_name,
  ...,
  ingest_mode = "write",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  soma_context = NULL
) {
  check_package("SummarizedExperiment", "1.28.0")
  stopifnot(
    "'uri' must be a single character value" = is.null(uri) ||
      is_scalar_character(uri),
    "'ms_name' must be a single character value" = is_scalar_character(
      ms_name
    ) &&
      nzchar(ms_name) &&
      !is.na(ms_name)
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

  experiment <- SOMAExperimentCreate(
    uri = uri,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    soma_context = soma_context
  )
  on.exit(experiment$close(), add = TRUE, after = FALSE)

  # Write cell-level meta data (obs)
  soma_info("Adding colData")
  obs_df <- .df_index(SummarizedExperiment::colData(x), axis = "obs")
  obs_df[[attr(obs_df, "index")]] <- colnames(x)
  write_soma(
    x = obs_df,
    uri = "obs",
    soma_parent = experiment,
    key = "obs",
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    soma_context = soma_context
  )

  # Write assays
  soma_info("Writing assays")
  expms <- SOMACollectionCreate(
    file_path(experiment$uri, "ms"),
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    soma_context = soma_context
  )
  withCallingHandlers(
    expr = .register_soma_object(expms, soma_parent = experiment, key = "ms"),
    existingKeyWarning = .maybe_muffle
  )
  ms_uri <- .check_soma_uri(uri = ms_name, soma_parent = expms)
  ms <- SOMAMeasurementCreate(
    uri = ms_uri,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    soma_context = soma_context
  )
  on.exit(ms$close(), add = TRUE, after = FALSE)

  X <- if (!"X" %in% ms$names()) {
    SOMACollectionCreate(
      uri = file_path(ms$uri, "X"),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      soma_context = soma_context
    )
  } else {
    SOMACollectionOpen(file_path(ms$uri, "X"), mode = "WRITE")
  }
  withCallingHandlers(
    .register_soma_object(X, soma_parent = ms, key = "X"),
    existingKeyWarning = .maybe_muffle
  )
  on.exit(X$close(), add = TRUE, after = FALSE)

  for (assay in SummarizedExperiment::assayNames(x)) {
    soma_info(sprintf("Adding %s assay", assay))
    write_soma(
      x = SummarizedExperiment::assay(x, assay),
      uri = assay,
      soma_parent = X,
      sparse = TRUE,
      transpose = TRUE,
      key = assay,
      ingest_mode = ingest_mode,
      shape = rev(shape),
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      soma_context = soma_context
    )
  }

  # Write feature-level meta data
  soma_info("Adding rowData")
  var_df <- .df_index(SummarizedExperiment::rowData(x), axis = "var")
  if (!is.null(rownames(x))) {
    var_df[[attr(var_df, "index")]] <- rownames(x)
  }
  write_soma(
    x = var_df,
    uri = "var",
    soma_parent = ms,
    key = "var",
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    soma_context = soma_context
  )

  withCallingHandlers(
    expr = .register_soma_object(ms, soma_parent = expms, key = ms_name),
    existingKeyWarning = .maybe_muffle
  )

  return(experiment$uri)
}
