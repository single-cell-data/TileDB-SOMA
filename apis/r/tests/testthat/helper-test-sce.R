pbmc3k_sce <- function() {
  stopifnot(requireNamespace("pbmc3k", quietly = TRUE))
  e <- new.env()
  on.exit(expr = rm(e), add = TRUE, after = FALSE)
  load(
    file = system.file(
      "data",
      "pbmc3k.rda",
      package = "pbmc3k",
      mustWork = TRUE
    ),
    envir = e
  )
  load(
    file = system.file(
      "data",
      "pbmc3k.sce.logcounts.rda",
      package = "pbmc3k",
      mustWork = TRUE
    ),
    envir = e
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = e$pbmc3k[, colnames(e$pbmc3k.sce.logcounts)],
      logcounts = e$pbmc3k.sce.logcounts
    ),
    mainExpName = "RNA"
  )
  hvg <- system.file("extdata", "pbmc3k-sce", "hvg.txt", package = "tiledbsoma")
  if (nzchar(hvg)) {
    S4Vectors::metadata(sce) <- list(hvg = readLines(hvg))
  }
  clusters <- system.file(
    "extdata",
    "pbmc3k-sce",
    "clusters.txt",
    package = "tiledbsoma"
  )
  if (nzchar(clusters)) {
    SingleCellExperiment::colLabels(sce) <- as.factor(readLines(clusters))
  }
  decode <- function(uri) {
    stopifnot(nzchar(uri), file.exists(uri))
    x <- as.matrix(read.csv(
      uri,
      header = TRUE,
      row.names = 1L,
      colClasses = "character"
    ))
    return(matrix(
      as.numeric(x),
      nrow = nrow(x),
      ncol = ncol(x),
      dimnames = dimnames(x)
    ))
  }
  pca <- system.file("extdata", "pbmc3k-sce", "pca.csv", package = "tiledbsoma")
  if (nzchar(pca)) {
    SingleCellExperiment::reducedDim(sce, "PCA") <- decode(pca)
  }
  umap <- system.file(
    "extdata",
    "pbmc3k-sce",
    "umap.csv",
    package = "tiledbsoma"
  )
  if (nzchar(umap)) {
    SingleCellExperiment::reducedDim(sce, "UMAP") <- decode(umap)
  }
  return(sce)
}
