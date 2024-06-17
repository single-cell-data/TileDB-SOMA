#' @param sr SOMA read pointer
#' @param array Underlying \code{\link{SOMASparseNDArray}}
#' @param axis Axis to iterate over in a blockwise manner
#' @param size The size of each blockwise chunk to generate
#' @param reindex_disable_on_axis Additional axes that will not be re-indexed;
#' the following values may be used as shorthands for common settings:
#' \itemize{
#'  \item \dQuote{\code{TRUE}}: disable re-indexing on all axes
#'  \item \dQuote{\code{NA}}: re-index only on \code{axis}, disable
#'   re-indexing on all others
#'  \item \dQuote{\code{FALSE}}: re-index on \emph{all} axes, do \strong{not}
#'   disable re-indexing
#' }
