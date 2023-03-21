test_that("Basic mechanics", {
  uri <- withr::local_tempdir("soma-ms")

  measurement <- SOMAMeasurement$new(uri, internal_use_only = "allowed_use")

  expect_false(measurement$exists())
  expect_error(measurement$var, "Group does not exist.")

  measurement$create()
  # TODO: Determine behavior for retrieving empty obs/ms
  # expect_null(experiment$var)
  # expect_null(experiment$X)

  # Add var
  expect_error(measurement$var, "No member named 'var' found")

  var <- create_and_populate_var(file.path(uri, "var"))

  measurement$var <- var
  expect_equal(measurement$length(), 1)
  expect_true(inherits(measurement$var, "SOMADataFrame"))

  # Add X collection
  expect_error(measurement$X, "No member named 'X' found")
  expect_error(measurement$X <- var, "X must be a 'SOMACollection'")

  X <- SOMACollection$new(file.path(uri, "X"), internal_use_only = "allowed_use")
  X$create()

  measurement$X <- X
  expect_true(inherits(measurement$X, "SOMACollection"))
  expect_equal(measurement$length(), 2)

  # Add X layer
  expect_equal(measurement$X$length(), 0)
  nda <- create_and_populate_sparse_nd_array(file.path(uri, "X", "RNA"))
  measurement$X$set(nda, name = "RNA")
  expect_equal(measurement$X$length(), 1)
})
