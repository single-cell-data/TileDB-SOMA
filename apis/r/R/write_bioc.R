#' @rdname write_soma_objects
#'
#' @method write_soma DataFrame
#' @export
#'
write_soma.DataFrame <- function(
  x,
  uri,
  soma_parent,
  df_index = NULL,
  index_column_names = 'soma_joinid',
  ...,
  ingest_mode = 'write',
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  relative = TRUE
) {
  # Check for compound non-atomic/factor types
  for (i in names(x)) {
    if (!(is.atomic(x[[i]]) || is.factor(x[[i]]))) {
      stop("All columns in DataFrames must be atomic or factors", call. = FALSE)
    }
  }
  index <- attr(x, which = 'index')
  x <- suppressWarnings(as.data.frame(x), classes = "deprecatedWarning")
  attr(x, which = 'index') <- index
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
    relative = relative
  ))
}

#' @rdname write_soma_objects
#'
#' @method write_soma Hits
#' @export
#'
write_soma.Hits <- function(
  x,
  uri,
  soma_parent,
  sparse = TRUE,
  type = NULL,
  transpose = FALSE,
  ...,
  ingest_mode = 'write',
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
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
    relative = relative
  ))
}

#' Write a \code{\link[SingleCellExperiment:SingleCellExperiment-class]{SingleCellExperiment}}
#' object to a SOMA
#'
#' @inheritParams write_soma
#' @inheritParams write_soma_objects
#' @param ms_name Name for resulting measurement; defaults to
#' \code{\link[SingleCellExperiment]{mainExpName}(x)}
#'
#' @inherit write_soma.SummarizedExperiment return sections
#'
#' @section Writing Reduced Dimensions:
#' Reduced dimensions are written out as
#' \link[tiledbsoma:SOMASparseNDArray]{sparse matrices} within the `obsm` group
#' of \code{\link[tiledbsoma:SOMAMeasurement]{measurement}} names `ms_name`
#'
#' @section Writing Column Pairs:
#' Column-wise relationship matrices are written out as
#' \link[tiledbsoma:SOMASparseNDArray]{sparse matrices} within the `obsp` group
#' of \code{\link[tiledbsoma:SOMAMeasurement]{measurement}} names `ms_name`
#'
#' @section Writing Row Pairs:
#' Row-wise relationship matrices are written out as
#' \link[tiledbsoma:SOMASparseNDArray]{sparse matrices} within the `varp` group
#' of \code{\link[tiledbsoma:SOMAMeasurement]{measurement}} names `ms_name`
#'
#' @method write_soma SingleCellExperiment
#' @export
#'
write_soma.SingleCellExperiment <- function(
  x,
  uri,
  ms_name = NULL,
  ...,
  ingest_mode = 'write',
  platform_config = NULL,
  tiledbsoma_ctx = NULL
) {
  check_package('SingleCellExperiment', version = .MINIMUM_SCE_VERSION())
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
  ms_name <- ms_name %||% SingleCellExperiment::mainExpName(x)

  uri <- NextMethod(
    'write_soma',
    x,
    uri = uri,
    ms_name = ms_name,
    ...,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  experiment <- SOMAExperimentOpen(
    uri = uri,
    mode = 'WRITE',
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  on.exit(expr = experiment$close(), add = TRUE, after = FALSE)

  ms <- experiment$ms$get(ms_name)
  on.exit(ms$close(), add = TRUE, after = FALSE)

  # Write reduced dimensions
  spdl::info("Adding reduced dimensions")
  if (!'obsm' %in% ms$names()) {
    ms$obsm <- SOMACollectionCreate(
      uri = file.path(ms$uri, 'obsm'),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  } else {
    ms$obsm$reopen("WRITE")
  }
  for (rd in SingleCellExperiment::reducedDimNames(x)) {
    spdl::info("Adding reduced dimension {}", rd)
    write_soma(
      x = SingleCellExperiment::reducedDim(x, rd),
      uri = rd,
      soma_parent = ms$obsm,
      sparse = TRUE,
      key = rd,
      ingest_mode = ingest_mode,
      shape = if (is.null(shape)) {
        NULL
      } else {
        c(shape[2L], ncol(SingleCellExperiment::reducedDim(x, rd)))
      },
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  }

  # Write nearest-neighbor graphs
  if (!'obsp' %in% ms$names()) {
    ms$obsp <- SOMACollectionCreate(
      uri = file.path(ms$uri, 'obsp'),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  } else {
    ms$obsp$reopen("WRITE")
  }
  for (cp in SingleCellExperiment::colPairNames(x)) {
    spdl::info("Adding colPair {}", cp)
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
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  }

  # Write coexpression networks
  if (!'varp' %in% ms$names()) {
    ms$varp <- SOMACollectionCreate(
      uri = file.path(ms$uri, 'varp'),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  } else {
    ms$varp$reopen("WRITE")
  }
  for (rp in SingleCellExperiment::rowPairNames(x)) {
    spdl::info("Adding rowPair {}", rp)
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
      tiledbsoma_ctx = tiledbsoma_ctx
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
#' @param ms_name Name for resulting measurement
#'
#' @inherit write_soma return
#'
#' @section Writing `colData`:
#' `colData` is written out as a
#' \link[tiledbsoma:SOMADataFrame]{data frame} called \dQuote{`obs`} at
#' the \code{\link[tiledbsoma:SOMAExperiment]{experiment}} level
#'
#' @section Writing Assay Matrices:
#' Each \link[SummarizedExperiment:assay]{assay matrix} is written out as a
#' \link[tiledbsoma:SOMASparseNDArray]{sparse matrix} within the `X` group of
#' \code{\link[tiledbsoma:SOMAMeasurement]{measurement}} names `ms_name`. Names
#' for assay matrices within `X` are taken from the
#' \link[SummarizedExperiment:assayNames]{assay names}. Assay matrices are
#' transposed (samples as rows) prior to writing
#'
#' @section Writing `rowData`:
#' `rowData` is written out as a
#' \link[tiledbsoma:SOMADataFrame]{data frame} called \dQuote{`var`} at
#' the \code{\link[tiledbsoma:SOMAMeasurement]{measurement}} level
#'
#' @method write_soma SummarizedExperiment
#' @export
#'
write_soma.SummarizedExperiment <- function(
  x,
  uri,
  ms_name,
  ...,
  ingest_mode = 'write',
  platform_config = NULL,
  tiledbsoma_ctx = NULL
) {
  check_package('SummarizedExperiment', '1.28.0')
  stopifnot(
    "'uri' must be a single character value" = is.null(uri) ||
      is_scalar_character(uri),
    "'ms_name' must be a single character value" = is_scalar_character(ms_name) &&
      nzchar(ms_name) &&
      !is.na(ms_name)
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

  experiment <- SOMAExperimentCreate(
    uri = uri,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  on.exit(experiment$close(), add = TRUE, after = FALSE)

  # Write cell-level meta data (obs)
  spdl::info("Adding colData")
  obs_df <- .df_index(SummarizedExperiment::colData(x), axis = 'obs')
  obs_df[[attr(obs_df, 'index')]] <- colnames(x)
  experiment$obs <- write_soma(
    x = obs_df,
    uri = 'obs',
    soma_parent = experiment,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )

  # Write assays
  spdl::info("Writing assays")
  experiment$add_new_collection(
    object = SOMACollectionCreate(
      file_path(experiment$uri, 'ms'),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    ),
    key = 'ms'
  )
  ms_uri <- .check_soma_uri(uri = ms_name, soma_parent = experiment$ms)
  ms <- SOMAMeasurementCreate(
    uri = ms_uri,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  on.exit(ms$close(), add = TRUE, after = FALSE)

  if (!'X' %in% ms$names()) {
    ms$X <- SOMACollectionCreate(
      uri = file.path(ms$uri, 'X'),
      ingest_mode = ingest_mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  } else {
    ms$X$reopen("WRITE")
  }
  on.exit(ms$X$close(), add = TRUE, after = FALSE)

  for (assay in SummarizedExperiment::assayNames(x)) {
    spdl::info("Adding {} assay", assay)
    write_soma(
      x = SummarizedExperiment::assay(x, assay),
      uri = assay,
      soma_parent = ms$X,
      sparse = TRUE,
      transpose = TRUE,
      key = assay,
      ingest_mode = ingest_mode,
      shape = rev(shape),
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
  }

  # Write feature-level meta data
  spdl::info("Adding rowData")
  var_df <- .df_index(SummarizedExperiment::rowData(x), axis = 'var')
  if (!is.null(rownames(x))) {
    var_df[[attr(var_df, 'index')]] <- rownames(x)
  }
  ms$var <- write_soma(
    x = var_df,
    uri = 'var',
    soma_parent = ms,
    ingest_mode = ingest_mode,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )

  experiment$ms$set(object = ms, name = ms_name)

  return(experiment$uri)
}
