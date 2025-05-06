## code to prepare `model-gene-var` dataset goes here

stopifnot(
  requireNamespace("scran", quietly = TRUE),
  requireNamespace("scater", quietly = TRUE),
  requireNamespace("bluster", quietly = TRUE),
  requireNamespace("pbmc3k", quietly = TRUE),
  requireNamespace("rprojroot", quietly = TRUE)
)

set.seed(42L)

extdata <- file.path(
  rprojroot::find_package_root_file(),
  "inst",
  "extdata",
  "pbmc3k-sce"
)
if (dir.exists(extdata)) {
  unlink(extdata, recursive = TRUE, force = TRUE)
}
dir.create(extdata, showWarnings = FALSE, recursive = TRUE)

data("pbmc3k.sce.logcounts", package = "pbmc3k")
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(logcounts = pbmc3k.sce.logcounts)
)

dec <- scran::modelGeneVar(sce)
hvg <- scran::getTopHVGs(dec, prop = 0.1)

sce <- scater::runPCA(sce, ncomponents = 50, subset_row = hvg)
sce <- scater::runUMAP(sce, dimred = "PCA")
clusters <- scran::clusterCells(
  sce,
  use.dimred = "PCA",
  BLUSPARAM = bluster::NNGraphParam(cluster.fun = "louvain")
)

writeLines(hvg, file.path(extdata, "hvg.txt"))
writeLines(as.vector(clusters), file.path(extdata, "clusters.txt"))

encode <- function(x) {
  stopifnot(is.matrix(x), is.double(x))
  return(matrix(
    data = sprintf(fmt = "%a", x),
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dimnames(x)
  ))
}

write.table(
  encode(SingleCellExperiment::reducedDim(sce, "PCA")),
  file = file.path(extdata, "pca.csv"),
  quote = FALSE,
  sep = ",",
  row.names = TRUE,
  col.names = TRUE
)
write.table(
  encode(SingleCellExperiment::reducedDim(sce, "UMAP")),
  file = file.path(extdata, "umap.csv"),
  quote = FALSE,
  sep = ",",
  row.names = TRUE,
  col.names = TRUE
)
