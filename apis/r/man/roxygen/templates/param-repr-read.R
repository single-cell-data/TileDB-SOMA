#' @param repr Representation of the sparse matrix to return; choose from:
#' \itemize{
#'  \item \dQuote{\code{T}}: returns a
#'   \code{\link[Matrix:TsparseMatrix-class]{TsparseMatrix}}
#'  \item \dQuote{\code{R}}: returns an
#'   \code{\link[Matrix:RsparseMatrix-class]{RsparseMatrix}}
#'  \item \dQuote{\code{C}}: returns a
#'   \code{\link[Matrix:CsparseMatrix-class]{CsparseMatrix}}
#' }
#' \strong{Note}: passing \code{repr} of \dQuote{\code{R}} or \dQuote{\code{C}}
#' are only available if re-indexing is enabled on axes \code{0} or \code{1},
#' respectively
