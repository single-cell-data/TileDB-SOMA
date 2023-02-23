is_seurat_assay <- function(x) {
  inherits(x, "Assay")
}

seurat_assay_has_scale_data <- function(x) {
  stopifnot(is_seurat_assay(x))
  !SeuratObject::IsMatrixEmpty(SeuratObject::GetAssayData(x, 'scale.data'))
  # length(SeuratObject::GetAssayData(x, "scale.data")) > 0
}
