#' Seurat 3k PBMCs from 10x Genomics
#'
#' Create a [`SeuratObject::Seurat`] object containing the widely used 3k PBMCs
#' dataset from 10x Genomics.
#'
#' @returns a [`SeuratObject::Seurat`] object
#' @seealso https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k
#' @importFrom Matrix readMM
#' @importFrom utils download.file read.table untar
#' @export

dataset_seurat_pbmc3k <- function() {
  url <- "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"

  tarfile <- file.path(tempdir(), basename(url))
  tardir <- sub("\\.tar\\.gz$", "", tarfile)
  datadir <- file.path(tardir, "filtered_gene_bc_matrices", "hg19")

  if (!dir.exists(tardir)) {
    if (!file.exists(tarfile)) {
      utils::download.file(url = url, destfile = tarfile)
    }
    utils::untar(tarfile, exdir = tardir)
  }

  mat <- Matrix::readMM(file.path(datadir, "matrix.mtx"))
  genes <- utils::read.table(
    file = file.path(datadir, "genes.tsv"),
    header = FALSE,
    col.names = c("id", "gene_name")
  )
  barcodes <- utils::read.table(
    file = file.path(datadir, "barcodes.tsv"),
    header = FALSE,
    col.names = "id"
  )
  dimnames(mat) <- list(genes$id, barcodes$id)

  object <- SeuratObject::CreateSeuratObject(counts = mat)
  object[["RNA"]] <- SeuratObject::AddMetaData(
    object = object[["RNA"]],
    metadata = genes$gene_name,
    col.name = "gene_name"
  )
  object
}
