# Create SOMAExperiment for the pbmc_small dataset
# https://mojaveazure.github.io/seurat-object/reference/pbmc_small.html

library(tiledbsoma)
library(SeuratObject)

# load the pbmc_small dataset
data(pbmc_small, package = "SeuratObject")
pbmc_small

# variables
data_dir <- normalizePath(file.path("inst", "extdata"))
soma_exp_name <- "soma-exp-pbmc-small"
soma_exp_uri <- file.path(tempdir(), soma_exp_name)
tar_file <- file.path(data_dir, paste0(soma_exp_name, ".tar.gz"))

# create the SOMAExperiment
write_soma(pbmc_small, uri = soma_exp_uri)

# create tar.gz file containing the SOMAExperiment
od <- setwd(soma_exp_uri)
tar(tar_file, compression = "gzip")
setwd(od)
