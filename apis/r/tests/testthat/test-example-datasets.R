test_that("example dataset access", {
  expect_length(
    list_datasets(),
    length(dir(example_data_dir()))
  )

  # Test that the dataset can be extracted
  dataset_uri <- extract_dataset("soma-exp-pbmc-small")
  expect_true(dir.exists(dataset_uri))
  expect_equal(tiledb::tiledb_object_type(dataset_uri), "GROUP")

  # Test the datasets can be loaded
  exp <- load_dataset("soma-exp-pbmc-small")
  exp <- SOMAExperimentOpen(exp$uri)
  expect_s3_class(exp, "SOMAExperiment")
  exp$close()

  sdf <- load_dataset("soma-dataframe-pbmc3k-processed-obs")
  sdf <- SOMADataFrameOpen(sdf$uri)
  expect_s3_class(sdf, "SOMADataFrame")
  sdf$close()
})
