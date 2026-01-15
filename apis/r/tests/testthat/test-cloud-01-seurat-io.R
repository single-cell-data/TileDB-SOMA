# Tests for TileDB Cloud Integration with Seurat I/O ----------------------

test_that("Ingest Seurat object to cloud with write_soma", {
  skip_if_no_cloud()
  skip_if_not_installed("SeuratObject", .MINIMUM_SEURAT_VERSION("c"))

  # Get test Seurat object
  pbmc_small <- get_test_seurat_object()

  # Create cloud experiment path
  uri <- cloud_path()

  # Ingest Seurat object to cloud --------------------------------------------
  expect_no_error(result_uri <- write_soma(pbmc_small, uri = uri))
  expect_equal(result_uri, uri)

  # Verify experiment exists and has correct metadata
  exp <- SOMAExperimentOpen(uri)
  expect_true(exp$exists())
  expect_equivalent(exp$get_metadata("dataset_type"), "soma")
  exp$close()

  # Access SOMA components from cloud experiment ----------------------------
  exp <- SOMAExperimentOpen(uri)

  # Verify obs dataframe
  expect_true(exp$obs$exists())
  expect_equal(exp$obs$soma_type, "SOMADataFrame")
  obs_data <- exp$obs$read()$concat()
  expect_s3_class(obs_data, "Table")
  expect_true("soma_joinid" %in% names(obs_data))
  expect_true("nFeature_RNA" %in% names(obs_data))

  # Verify measurement
  ms_rna <- exp$ms$get("RNA")
  expect_equivalent(ms_rna$soma_type, "SOMAMeasurement")

  # Verify var dataframe
  var_data <- ms_rna$var$read()$concat()
  expect_s3_class(var_data, "Table")
  expect_true("soma_joinid" %in% names(var_data))
  expect_true("var_id" %in% names(var_data))

  # Slice and filter obs dataframe with coords ------------------------------
  exp <- SOMAExperimentOpen(uri)

  # Test filtering with coords
  obs_coords <- exp$obs$read(
    coords = 0:9,
    column_names = c("soma_joinid", "obs_id", "nFeature_RNA")
  )$concat()
  expect_equal(nrow(obs_coords), 10)

  # Test filtering with value_filter
  obs_filtered <- exp$obs$read(
    value_filter = "nFeature_RNA > 70",
    column_names = c("soma_joinid", "obs_id", "nFeature_RNA")
  )$concat()$to_data_frame()
  expect_true(nrow(obs_filtered) > 0)
  # Verify all rows meet the filter criteria
  expect_true(all(obs_filtered$nFeature_RNA > 70))

  # Test combining coords and value_filter
  obs_combined <- exp$obs$read(
    coords = 0:49,
    value_filter = "nFeature_RNA > 70",
    column_names = c("soma_joinid", "obs_id", "nFeature_RNA")
  )$concat()$to_data_frame()

  expect_true(nrow(obs_combined) <= 50)
  expect_true(all(obs_combined$nFeature_RNA > 70))
  exp$close()

  # Slice and filter var dataframe ------------------------------------------
  exp <- SOMAExperimentOpen(uri)
  ms_rna <- exp$ms$get("RNA")

  # Get a few gene names from the var dataframe
  all_var <- ms_rna$var$read()$concat()$to_data_frame()
  test_genes <- head(all_var$var_id, 4)

  # Test filtering with gene names using IN operator
  var_filter_str <- sprintf(
    "var_id %%in%% c(%s)",
    paste0("'", test_genes, "'", collapse = ", ")
  )
  var_filtered <- ms_rna$var$read(
    value_filter = var_filter_str,
    column_names = c("soma_joinid", "var_id")
  )$concat()$to_data_frame()

  expect_equal(nrow(var_filtered), length(test_genes))
  expect_setequal(var_filtered$var_id, test_genes)
  exp$close()

  # SOMAExperimentAxisQuery for filtering experiments
  exp <- SOMAExperimentOpen(uri)

  # Create query with obs and var filtering
  query <- SOMAExperimentAxisQuery$new(
    experiment = exp,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(
      value_filter = "nFeature_RNA > 70"
    ),
    var_query = SOMAAxisQuery$new(
      coords = 0:99
    )
  )

  # Verify query dimensions
  expect_true(query$n_obs > 0)
  expect_equal(query$n_vars, 100)

  # Read obs from query
  obs_result <- query$obs(column_names = c("obs_id", "nFeature_RNA"))$concat()
  expect_equal(nrow(obs_result), query$n_obs)
  expect_true(all(obs_result$nFeature_RNA > 70))

  # Read X data from query
  x_reader <- query$X(layer_name = "data")
  expect_s3_class(x_reader, "SOMASparseNDArrayRead")

  # Convert to Seurat
  seurat_result <- query$to_seurat(X_layers = c(data = "data"))
  expect_s3_class(seurat_result, "Seurat")
  expect_equal(ncol(seurat_result), query$n_obs)
  expect_equal(nrow(seurat_result), query$n_vars)
  exp$close()

  # Add dimensional reduction to cloud experiment ---------------------------
  exp <- SOMAExperimentOpen(uri, mode = "READ")
  ms_rna <- exp$ms$get("RNA")

  # obsm may not exist or be empty when no reductions are present
  if (ms_rna$obsm$exists()) {
    expect_false("PCA" %in% ms_rna$obsm$names())
  }
  exp$close()

  # Now prepare the object with PCA
  suppressWarnings({
    pbmc_small <- Seurat::FindVariableFeatures(
      pbmc_small,
      nfeatures = 50,
      verbose = FALSE
    )
    pbmc_small <- Seurat::ScaleData(pbmc_small, verbose = FALSE)
    pbmc_small <- Seurat::RunPCA(
      pbmc_small,
      features = Seurat::VariableFeatures(pbmc_small),
      verbose = FALSE
    )
  })

  # Reopen experiment in WRITE mode
  exp <- SOMAExperimentOpen(uri, mode = "WRITE")
  expect_equal(exp$mode(), "WRITE")

  # Get RNA measurement in WRITE mode
  ms_rna <- exp$ms$get("RNA")
  expect_equal(ms_rna$mode(), "WRITE")

  # Add PCA to the experiment
  suppressWarnings({
    write_soma(pbmc_small[["pca"]], soma_parent = ms_rna)
  })

  # Close measurement and experiment
  ms_rna$close()
  exp$close()

  # Reopen in READ mode and verify PCA exists
  exp <- SOMAExperimentOpen(uri, mode = "READ")
  ms_rna <- exp$ms$get("RNA")

  # Verify obsm collection exists with PCA
  expect_true(ms_rna$obsm$exists())
  expect_true("PCA" %in% ms_rna$obsm$names())

  # Verify PCA data can be read
  pca_data <- ms_rna$obsm$get("PCA")
  expect_equivalent(pca_data$soma_type, "SOMADenseNDArray")
  pca_array <- pca_data$read()$concat()
  expect_s3_class(pca_array, "Table")
  exp$close()

  # Filter obs with string/dictionary column equality -----------------------
  exp <- SOMAExperimentOpen(uri)

  # Get total row count
  total_obs <- exp$obs$read()$concat()
  total_rows <- nrow(total_obs)

  # Test string equality filter on orig.ident (dictionary column)
  # All cells in pbmc_small have orig.ident == "SeuratProject"
  obs_equal <- exp$obs$read(
    coords = 0:9,
    value_filter = 'orig.ident == "SeuratProject"',
    column_names = c("soma_joinid", "obs_id", "orig.ident")
  )$concat()
  expect_equal(nrow(obs_equal), 10)

  # Test string inequality filter
  # (should return 0 rows since all are SeuratProject)
  obs_not_equal <- exp$obs$read(
    coords = 0:9,
    value_filter = 'orig.ident != "SeuratProject"',
    column_names = c("soma_joinid", "obs_id", "orig.ident")
  )$concat()
  expect_equal(nrow(obs_not_equal), 0)
  exp$close()

  # Data integrity after round-trip to cloud --------------------------------
  # Capture original dimensions and cell IDs
  original_ncells <- ncol(pbmc_small)
  original_ngenes <- nrow(pbmc_small)
  original_cell_ids <- colnames(pbmc_small)
  original_gene_ids <- rownames(pbmc_small)

  # Write to cloud
  write_soma(pbmc_small, uri = uri)

  # Read back
  exp <- SOMAExperimentOpen(uri)

  # Verify obs dimensions
  obs_data <- exp$obs$read()$concat()$to_data_frame()
  expect_equal(nrow(obs_data), original_ncells)
  expect_setequal(obs_data$obs_id, original_cell_ids)

  # Verify var dimensions
  ms_rna <- exp$ms$get("RNA")
  var_data <- ms_rna$var$read()$concat()$to_data_frame()
  expect_equal(nrow(var_data), original_ngenes)
  expect_setequal(var_data$var_id, original_gene_ids)

  # Verify X layer dimensions via sparse read
  x_data <- ms_rna$X$get("data")
  x_shape <- x_data$shape()
  expect_equal(x_shape[1], original_ncells)
  expect_equal(x_shape[2], original_ngenes)

  exp$close()
})

test_that("Error handling for cloud operations", {
  skip_if_no_cloud()

  # Test 1: Opening non-existent URI should error
  nonexistent_uri <- file_path(get_cloud_base_uri(), "does_not_exist_12345")
  expect_error(
    SOMAExperimentOpen(nonexistent_uri),
    class = "error"
  )

  # Test 2: Operations on closed objects should error
  skip_if_not_installed("SeuratObject", .MINIMUM_SEURAT_VERSION("c"))

  pbmc_small <- get_test_seurat_object()
  uri <- cloud_path()
  write_soma(pbmc_small, uri = uri)

  exp <- SOMAExperimentOpen(uri, mode = "READ")
  exp$close()

  # Attempting to use closed experiment should error
  expect_error(
    exp$obs,
    class = "error"
  )
})

test_that("SOMACollection member operations in cloud", {
  skip_if_no_cloud()
  uri <- cloud_path()

  # Create collection with multiple member types
  SOMACollectionCreate(uri = uri)$close()

  collection <- SOMACollectionOpen(uri, mode = "WRITE")

  # Add different SOMA types as members
  schema <- arrow::schema(
    arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
    arrow::field("value", arrow::int32())
  )
  domain <- list(soma_joinid = c(0L, 100L))
  df <- collection$add_new_dataframe(
    "test_dataframe",
    schema = schema,
    index_column_names = "soma_joinid",
    domain = domain
  )
  df$close()

  sparse <- collection$add_new_sparse_ndarray(
    "test_sparse",
    type = arrow::float32(),
    shape = c(10L, 10L)
  )
  sparse$close()

  # Verify members exist
  expect_equal(collection$length(), 2)
  expect_setequal(
    collection$names(),
    c("test_dataframe", "test_sparse")
  )

  collection$close()

  # Verify members can be opened and have correct types
  collection <- SOMACollectionOpen(uri, mode = "READ")
  expect_equal(collection$get("test_dataframe")$soma_type, "SOMADataFrame")
  expect_equal(collection$get("test_sparse")$soma_type, "SOMASparseNDArray")
  collection$close()
})

test_that("DELETE mode member removal in cloud", {
  skip_if_no_cloud()

  uri <- cloud_path()

  # Create collection with members
  SOMACollectionCreate(uri = uri)$close()

  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  collection$add_new_sparse_ndarray(
    "array_to_delete",
    type = arrow::int8(),
    shape = c(10L)
  )$close()
  collection$add_new_sparse_ndarray(
    "array_to_keep",
    type = arrow::int8(),
    shape = c(10L)
  )$close()
  collection$close()

  # Verify both members exist
  collection <- SOMACollectionOpen(uri, mode = "READ")
  expect_setequal(collection$names(), c("array_to_delete", "array_to_keep"))
  collection$close()

  # Delete one member
  collection <- SOMACollectionOpen(uri, mode = "DELETE")
  expect_no_error(collection$remove("array_to_delete"))
  collection$close()

  # Verify only one member remains
  collection <- SOMACollectionOpen(uri, mode = "READ")
  expect_equal(collection$names(), "array_to_keep")
  expect_equal(collection$length(), 1)
  collection$close()

  # Test error when trying to remove non-existent member
  collection <- SOMACollectionOpen(uri, mode = "DELETE")
  expect_error(
    collection$remove("nonexistent"),
    "does not exist"
  )
  collection$close()
})
