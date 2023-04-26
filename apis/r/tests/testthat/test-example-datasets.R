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
  expect_s3_class(exp, "SOMAExperiment")

  sdf <- load_dataset("soma-dataframe-pbmc3k-processed-obs")
  expect_s3_class(sdf, "SOMADataFrame")
})
