# Create SOMADataFrame for pbmc3k_processed

library(tiledb)
library(tiledbsoma)
library(SeuratObject)

# Prerequisite:
#
#   install.packages("pbmc3k.SeuratData", repos = c("https://seurat.nygenome.org", getOption("repos")))
#
# Note: as of 2024-11-07:
#
#   Warning: unable to access index for repository https://cran.r-project.org/src/contrib:
#     download from 'https://cran.r-project.org/src/contrib/PACKAGES' failed
#   Warning: unable to access index for repository https://seurat.nygenome.org/bin/macosx/big-sur-arm64/contrib/4.3:
#     cannot open URL 'https://seurat.nygenome.org/bin/macosx/big-sur-arm64/contrib/4.3/PACKAGES'
#
# Advice: run this on a Linux system
library(pbmc3k.SeuratData)

data(pbmc3k.final)
pbmc3k.final <- UpdateSeuratObject(pbmc3k.final)

data_dir <- normalizePath(file.path("inst", "extdata"))

soma_exp_name <- "soma-dataframe-pbmc3k-processed"
soma_exp_uri  <- file.path(tempdir(), soma_exp_name)

# Create the SOMAExperiment
write_soma(pbmc3k.final, uri = soma_exp_uri)

soma_obs_name <- "soma-dataframe-pbmc3k-processed-obs"
soma_obs_uri <- paste(soma_exp_uri, "obs", sep="/")
tar_file <- file.path(data_dir, paste0(soma_obs_name, ".tar.gz"))

# Create tar.gz file containing the SOMAExperiment
od <- setwd(soma_obs_uri)
tar(tar_file, compression = "gzip")
setwd(od)
