library(tiledbsc)
library(lobstr)
library(tiledb)
library(SeuratObject)

# devtools::install_github('satijalab/seurat-data')
library(SeuratData)

# parameters -------------------------------------------------------------------

# where the new arrays will be stored
array_dir <- "dev/data/arrays"

# where the benchmark results will be stored
output_dir <- "dev/data/benchmarks"

# should the new arrays be deleted after the benchmark?
clean <- FALSE

# should datasets be installed if they are not already?
install_missing <- TRUE


# outputs ----------------------------------------------------------------------

# create the output file
output_metadata <- sprintf(
  "%s_tiledbsc%s_tiledb%s",
  format(Sys.time(),"%Y%m%d-%H%M"),
  packageVersion("tiledbsc"),
  packageVersion("tiledb")
)

output_name <- paste0("ingestion-benchmarks_", output_metadata, ".csv")

output_file <- file.path(output_dir, output_name)

# create parent directory using file metadata
array_dir <- file.path(array_dir, output_metadata)
dir.create(array_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# setup datasets ---------------------------------------------------------------
datasets <- AvailableData()

# remove datasets with Versions < 2.0.0 to avoid the following error:
# Not a validObject(): no slot of name "images" for this object of class "Seurat"
datasets <- subset(datasets, Version >= package_version("2.0.0"))

# skip datasets that cause errors
blocklist <- c(
  # Not a validObject(): no slot of name "images" for this object of class "Seurat"
  # I believe this error is specific to datasets with Versions < 2.0.0
  "humancortexref.SeuratData",
  "kidneyref.SeuratData",

  # Datasets that are already installed but LoadData() errors with:
  # Error: Could not find dataset '<>', please check manifest and try again
  "bonemarrowref.SeuratData",

  # Failed to download
  "lungref.SeuratData"
)

# Temporarily subset
datasets <- datasets[1:min(10L, nrow(datasets)),]

# setup results ----------------------------------------------------------------

# initialize results dataframe
benchmarks <- data.frame(
  size_memory = numeric(),
  size_rds = numeric(),
  size_tiledb = numeric(),
  ingest_time = numeric()
)

# main -------------------------------------------------------------------------
for (ds_name in rownames(datasets)) {
  message(sprintf("Dataset: %s", ds_name))
  if (ds_name %in% blocklist) {
    message(sprintf("..skipping blocklisted dataset: %s", ds_name))
    next
  }

  if (!datasets[ds_name, "Installed"]) {
    skip <- FALSE
    if (interactive()) {
      if (menu(c("yes", "no"), sprintf("Install %s?", ds_name)) == 2) skip <- TRUE
    } else {
      if (!install_missing) skip <- TRUE
    }

    if (skip) {
      message(sprintf("..skipping dataset not installed: %s", ds_name))
      next
    }

    message(sprintf("..installing dataset '%s'", ds_name))
    install_worked <- try(InstallData(ds_name), silent = TRUE)
    if (inherits(install_worked, "try-error")) {
      message(sprintf("..failed to install dataset '%s'", ds_name))
      next
    }
  }

  message("..loading dataset")
  ds <- LoadData(ds_name)

  uri <- file.path(array_dir, ds_name)
  message(sprintf("..ingesting data into '%s'", uri))

  ingest_start <- Sys.time()
  scgroup <- SCGroup$new(uri = uri)
  scgroup$from_seurat_assay(
    object = ds[[DefaultAssay(ds)]],
    obs = ds[[]]
  )

  benchmarks[ds_name, ] <- list(
    size_memory = as.numeric(obj_size(ds)),
    size_rds = tiledb_vfs_dir_size(system.file("data", package = ds_name)),
    size_tiledb = tiledb_vfs_dir_size(uri),
    ingest_time = as.numeric(difftime(Sys.time(), ingest_start, units = "secs"))
  )

  # write each time to avoid losing data on failures
  write.csv(benchmarks, file = output_file, quote = FALSE)

  if (clean) {
    message(sprintf("..cleaning dataset '%s'", ds_name))
    tiledb_vfs_remove_dir(uri)
  }
}
