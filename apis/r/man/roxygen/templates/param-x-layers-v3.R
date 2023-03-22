#' @param X_layers A named character of X layers to add to the Seurat assay,
#' where the names are the names of Seurat slots and the values are the names
#' of layers within \code{X}; names should be one of:
#' \itemize{
#'  \item \dQuote{\code{counts}} to add the layer as \code{counts}
#'  \item \dQuote{\code{data}} to add the layer as \code{data}
#'  \item \dQuote{\code{scale.data}} to add the layer as \code{scale.data}
#' }
#' At least one of \dQuote{\code{counts}} or \dQuote{\code{data}} is required
