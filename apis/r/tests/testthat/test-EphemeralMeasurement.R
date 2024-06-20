test_that("Ephemeral Measurement mechanics", {
  # Create the measurement
  uri <- tempfile(pattern="ephemeral-ms")
  expect_warning(EphemeralMeasurement$new(uri))
  expect_no_condition(measurement <- EphemeralMeasurement$new())
  expect_true(grepl('^ephemeral-collection:0x[[:digit:]a-f]{6,32}$', measurement$uri))
  expect_false(measurement$exists())
  expect_error(measurement$var)
  expect_s3_class(measurement$create(), measurement$class())

  # Add var
  expect_error(measurement$var)
  var <- create_and_populate_var(file.path(uri, "var"))
  measurement$var <- var
  expect_equal(measurement$length(), 1)
  expect_s3_class(measurement$var, 'SOMADataFrame')

  # Add X collection
  expect_error(measurement$X)
  expect_error(measurement$X <- var)
  X <- SOMACollectionCreate(file.path(uri, "X"))
  expect_no_condition(measurement$X <- X)
  expect_s3_class(measurement$X, 'SOMACollection')
  expect_equal(measurement$length(), 2)

  # Add X layer
  expect_equal(measurement$X$length(), 0)
  nda <- create_and_populate_sparse_nd_array(file.path(uri, "X", "RNA"))
  expect_no_condition(measurement$X$set(nda, name = "RNA"))
  expect_equal(measurement$X$length(), 1)

  X$close()
})
