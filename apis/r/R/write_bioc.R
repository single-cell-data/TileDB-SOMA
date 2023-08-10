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
  x <- as.data.frame(x)
  attr(x, which = 'index') <- index
  return(write_soma(
    x = x,
    uri = uri,
    soma_parent = soma_parent,
    df_index = df_index,
    index_column_names = index_column_names,
    ...,
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
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    relative = relative
  ))
}

#' Write a \code{\link[SingleCellExperiment:SingleCellExperiment-class]{SingleCellExperiment}}
#' object to a SOMA
#'
#' @inheritParams write_soma
#' @param ms_name Name for resulting measurement; defaults to
#' \code{\link[SingleCellExperiment]{mainExpName}(x)}
#'
#' @inherit write_soma return
#'
#' @inherit write_soma.SummarizedExperiment sections
#'
#' @section Writing Reduced Dimensions:
#' blah
#'
#' @section writing Column Pairs:
#' blah
#'
#' @section Writing Row Pairs:
#' blah
#'
#' @method write_soma SingleCellExperiment
#' @export
#'
write_soma.SingleCellExperiment <- function(
  x,
  uri,
  ms_name = NULL,
  ...,
  platform_config = NULL,
  tiledbsoma_ctx = NULL
) {
  check_package('SingleCellExperiment', version = .MINIMUM_SCE_VERSION())
  ms_name <- ms_name %||% SingleCellExperiment::mainExpName(x)
  uri <- NextMethod(
    'write_soma',
    x,
    uri = uri,
    ms_name = ms_name,
    ...,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  experiment <- SOMAExperimentOpen(
    uri = uri,
    mode = 'WRITE',
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  on.exit(expr = experiment$close(), add = TRUE)
  ms <- experiment$ms$get(ms_name)

    # Write reduced dimensions
  spdl::info("Adding reduced dimensions")
  obsm <- tryCatch(
    expr = ms$obsm,
    error = function(...) {
      return(NULL)
    }
  )
  if (is.null(obsm)) {
    ms$obsm <- SOMACollectionCreate(
      uri = file.path(ms$uri, 'obsm'),
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
    obsm <- ms$obsm
  }
  for (rd in SingleCellExperiment::reducedDimNames(x)) {
    spdl::info("Adding reduced dimension {}", rd)
    obsm$set(
      object = write_soma(
        x = SingleCellExperiment::reducedDim(x, rd),
        uri = rd,
        soma_parent = obsm,
        sparse = TRUE,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      ),
      name = rd
    )
  }

  # Write nearest-neighbor graphs
  obsp <- tryCatch(
    expr = ms$obsp,
    error = function(...) {
      return(NULL)
    }
  )
  if (is.null(obsp)) {
    ms$obsp <- SOMACollectionCreate(
      uri = file.path(ms$uri, 'obsp'),
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
    obsp <- ms$obsp
  }
  for (cp in SingleCellExperiment::colPairNames(x)) {
    spdl::info("Adding colPair {}", cp)
    obsp$set(
      object = write_soma(
        x = SingleCellExperiment::colPair(x, cp),
        uri = cp,
        soma_parent = obsp,
        sparse = TRUE,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      ),
      name = cp
    )
  }

  # Write coexpression networks
  varp <- tryCatch(
    expr = ms$varp,
    error = function(...) {
      return(NULL)
    }
  )
  if (is.null(varp)) {
    ms$varp <- SOMACollectionCreate(
      uri = file.path(ms$uri, 'varp'),
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    )
    varp <- ms$varp
  }
  for (rp in SingleCellExperiment::rowPairNames(x)) {
    spdl::info("Adding rowPair {}", rp)
    varp$set(
      object = write_soma(
        x = SingleCellExperiment::rowPair(x, rp),
        uri = rp,
        soma_parent = varp,
        sparse = TRUE,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      ),
      name = rp
    )
  }

  # TODO: Add alternate experiments
  return(experiment$uri)
}

#' Write a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object to a SOMA
#'
#' @inheritParams write_soma
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
  experiment <- SOMAExperimentCreate(
    uri = uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  on.exit(experiment$close(), add = TRUE)

  # Write cell-level meta data (obs)
  spdl::info("Adding colData")
  obs_df <- .df_index(SummarizedExperiment::colData(x), axis = 'obs')
  obs_df[[attr(obs_df, 'index')]] <- colnames(x)
  experiment$obs <- write_soma(
    x = obs_df,
    uri = 'obs',
    soma_parent = experiment,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )

  # Write assays
  spdl::info("Writing assays")
  experiment$add_new_collection(
    object = SOMACollectionCreate(
      file_path(experiment$uri, 'ms'),
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx
    ),
    key = 'ms'
  )
  ms_uri <- .check_soma_uri(uri = ms_name, soma_parent = experiment$ms)
  ms <- SOMAMeasurementCreate(
    uri = ms_uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  ms$X <- SOMACollectionCreate(
    uri = file.path(ms$uri, 'X'),
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  for (assay in SummarizedExperiment::assayNames(x)) {
    spdl::info("Adding {} assay", assay)
    ms$X$set(
      object = write_soma(
        x = SummarizedExperiment::assay(x, assay),
        uri = assay,
        soma_parent = ms$X,
        sparse = TRUE,
        transpose = TRUE,
        platform_config = platform_config,
        tiledbsoma_ctx = tiledbsoma_ctx
      ),
      name = assay
    )
  }
  ms$X$close()

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
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )

  ms$close()
  experiment$ms$set(object = ms, name = ms_name)

  return(experiment$uri)
}
