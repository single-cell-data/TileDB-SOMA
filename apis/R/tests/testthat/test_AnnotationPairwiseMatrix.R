test_that("annotation pairwise matrix can be stored and retrieved", {

  uri <- file.path(withr::local_tempdir(), "annot-mat")
  mat <- matrix(runif(100), nrow = 10, ncol = 10)

  annotmat <- AnnotationPairwiseMatrix$new(uri)
  expect_true(inherits(annotmat, "AnnotationPairwiseMatrix"))

  rlabs <- clabs <- paste0("R", seq_len(nrow(mat)))
  dimnames(mat) <- list(rlabs, clabs)

  annotmat$from_matrix(mat, index_cols = c("i", "j"))
  expect_true(dir.exists(annotmat$uri))
  expect_s4_class(annotmat$tiledb_array(), "tiledb_array")

  mat2 <- annotmat$to_matrix()
  expect_identical(mat[rlabs, clabs], mat2[rlabs, clabs])
})
