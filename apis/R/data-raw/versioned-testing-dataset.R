# Create an SCDataset containing pbmc_small using the current version of
# tiledbsc for backwards compatibility testing

library(tiledb)
library(tiledbsc)
library(SeuratObject)
data("pbmc_small", package = "SeuratObject")

data_dir <- "tests/testthat/test-data"

# version specific dataset name
dataset_name <- paste(
  "scdataset-pbmc-small",
  gsub("\\.", "-", packageVersion("tiledbsc")),
  sep = "_"
)

# setup directory
uri <- file.path(data_dir, dataset_name)

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
if (dir.exists(uri)) unlink(uri, recursive = TRUE)

# ingest
scdataset <- SCDataset$new(uri, verbose = TRUE)
scdataset$from_seurat(pbmc_small)

# populate the SCGroup's misc slot with an additional copy of obs so
# it's not empty and we cover the use case of a manually added misc item
scgroup <- scdataset$members$RNA
obs2_uri <- file.path(scgroup$misc$uri, "RNA_obs")
fromDataFrame(obj = pbmc_small[[]], uri = obs2_uri)
obs2 <- TileDBArray$new(obs2_uri)
scgroup$misc$add_member(obs2, name = "RNA_obs", relative = TRUE)

# zip the entire array to avoid build problems on windows caused by path length
# limitations we use zip because ?tar warns that known problems arise from
# handling of file paths of more than 100 bytes
out <- zip(
  zipfile = paste0(uri, ".zip"),
  files = uri
)
