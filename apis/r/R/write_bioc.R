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
  .NotYetImplemented()
}

#' @method write_soma SingleCellExperiment
#' @export
#'
write_soma.SingleCellExperiment <- function(
  x,
  uri,
  ...,
  platform_config = NULL,
  tiledbsoma_ctx = NULL
) {
  uri <- NextMethod()
  experiment <- SOMAExperimentOpen(
    uri = uri,
    mode = 'WRITE',
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx
  )
  on.exit(expr = experiment$close(), add = TRUE)
  .NotYetImplemented()
  return(experiment$uri)
}

#' @method write_soma SummarizedExperiment
#' @export
#'
write_soma.SummarizedExperiment <- function(
  x,
  uri,
  ...,
  platform_config = NULL,
  tiledbsoma_ctx = NULL
) {
  .NotYetImplemented()
}
