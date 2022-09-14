test_that("New matrices can be added", {
  uri <- withr::local_tempdir("assay-matrix-group")
  imat1 <- create_sparse_matrix_with_string_dims()

  amats <- AssayMatrixGroup$new(uri = uri, dimension_name = c("i", "j"))
  expect_length(amats$members, 0)

  amats$add_assay_matrix(data = imat1, name = "mat1")

  amat1 <- amats$get_member("mat1")
  expect_equal(amat1$dimnames(), amats$dimension_name)

  # create an assaymat outside of the group with non-matching dimnames
  # amat2 <- AssayMatrix$new(file.path(amats$uri, "mat2"))
  # amat2$from_matrix(imat1, index_cols = c("I", "J"))

  # TODO - this should be an error
  # expect_error(
  #   amats$add_member(amat2, name = "mat2", relative = TRUE)
  # )
})

test_that("Transposed matrices can be added and retrieved", {
  uri <- withr::local_tempdir("transposed-assay-matrix-group")
  imat1 <- create_sparse_matrix_with_string_dims(seed = 1)
  imat2 <- create_sparse_matrix_with_string_dims(seed = 2)

  amats <- AssayMatrixGroup$new(uri = uri, dimension_name = c("j", "i"))
  amats$add_assay_matrix(data = imat1, name = "mat1", transpose = TRUE)

  amat1 <- amats$get_member("mat1")
  expect_equal(amat1$dimnames(), amats$dimension_name)

  # by default, the transpose is not applied on read
  omat1 <- amat1$to_matrix()
  expect_equal(dim(Matrix::t(omat1)), dim(imat1))

  # must set transpose=TRUE to restore the original matrix shape
  omat1 <- amat1$to_matrix(transpose = TRUE)
  expect_equal(dim(omat1), dim(imat1))
  expect_setequal(rownames(omat1), rownames(imat1))
  expect_setequal(colnames(omat1), colnames(imat1))
})
