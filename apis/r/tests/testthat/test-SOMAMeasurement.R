test_that("Basic mechanics", {
  uri <- withr::local_tempdir("soma-ms")

  measurement <- SOMAMeasurementCreate(uri)

  # TODO: Determine behavior for retrieving empty obs/ms
  # expect_null(experiment$ms)
  # expect_null(experiment$obs)

  # Add var
  expect_error(measurement$var, "No member named 'var' found")

  var <- create_and_populate_var(file.path(uri, "var"))

  measurement$var <- var
  expect_equal(measurement$length(), 1)
  expect_true(inherits(measurement$var, "SOMADataFrame"))

  # Add X collection
  expect_error(measurement$X, "No member named 'X' found")
  expect_error(measurement$X <- var, "X must be a 'SOMACollection'")

  measurement$X <- SOMACollectionCreate(file.path(uri, "X"))
  expect_equal(measurement$X$mode(), "WRITE")

  expect_equal(measurement$X$mode(), "WRITE")
  expect_true(inherits(measurement$X, "SOMACollection"))
  expect_equal(measurement$length(), 2)

  # Add X layer
  expect_equal(measurement$X$length(), 0)
  nda <- create_and_populate_sparse_nd_array(file.path(uri, "X", "RNA"))
  measurement$X$set(nda, name = "RNA")
  expect_equal(measurement$X$length(), 1)

  measurement$X$close()
  measurement$close()
})
