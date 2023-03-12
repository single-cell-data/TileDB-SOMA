test_that("annotation matrix can be stored and retrieved", {

  uri <- file.path(withr::local_tempdir(), "annot-mat")
  mat <- matrix(runif(1000), nrow = 100, ncol = 10)

  annotmat <- AnnotationMatrix$new(uri)
  expect_true(inherits(annotmat, "AnnotationMatrix"))
  expect_error(
    annotmat$from_matrix(mat, index_col = "obs_id"),
    "Annotation data must be a matrix with defined dim names"
  )

  rlabs <- paste0("R", seq_len(nrow(mat)))
  clabs <- paste0("C", seq_len(ncol(mat)))
  dimnames(mat) <- list(rlabs, clabs)

  annotmat$from_matrix(mat, index_col = "obs_id")
  expect_true(dir.exists(annotmat$uri))
  expect_s4_class(annotmat$tiledb_array(), "tiledb_array")

  mat2 <- annotmat$to_matrix()
  expect_equal(sort(rownames(mat2)), sort(rownames(mat)))
  expect_equal(sort(colnames(mat2)), sort(colnames(mat)))

  expect_identical(mat[rlabs, clabs], mat2[rlabs, clabs])

  # test that result is identical with batch mode
  mat3 <- annotmat$to_matrix(batch_mode = TRUE)
  expect_identical(mat2, mat3)

  # verify legacy metadata tag is present
  expect_equal(
    annotmat$get_metadata(SOMA_LEGACY_VALIDITY_KEY),
    SOMA_LEGACY_VALIDITY
  )
})
