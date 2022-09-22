# Load pbmc_small from SeuratObject
data("pbmc_small", package = "SeuratObject")


create_sparse_matrix_with_string_dims <- function(nrows = 10, ncols = 5, seed = 1, repr = "T") {
  set.seed(seed)
  smat <- Matrix::rsparsematrix(
    nrow = nrows,
    ncol = ncols,
    density = 0.6,
    rand.x = function(n) as.integer(runif(n, min = 1, max = 100)),
    repr = repr
  )
  dimnames(smat) <- list(
    paste0("i", seq_len(nrows)),
    paste0("j", seq_len(ncols))
  )
  smat
}
