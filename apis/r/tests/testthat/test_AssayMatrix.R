test_that("AssayMatrix object can be created from a dgCMatrix", {
  uri <- withr::local_tempdir("assay-matrix")
  mat <- SeuratObject::GetAssayData(pbmc_small[["RNA"]], "counts")

  assaymat <- AssayMatrix$new(uri = uri, verbose = FALSE)
  expect_true(inherits(assaymat, "AssayMatrix"))

  assaymat$from_matrix(mat, index_cols = c("i", "j"), value_col = "counts")
  expect_true(dir.exists(uri))
  expect_s4_class(assaymat$tiledb_array(), "tiledb_array")
  expect_equal(assaymat$dimnames(), c("i", "j"))

  df2 <- assaymat$to_dataframe()
  expect_s3_class(df2, "data.frame")
  expect_equal(attr(df2, "query_status"), "COMPLETE")
  expect_setequal(unique(df2$i), rownames(mat))
  expect_setequal(unique(df2$j), colnames(mat))

  mat2 <- assaymat$to_matrix()
  expect_s4_class(mat2, "dgTMatrix")
  expect_setequal(rownames(mat2), rownames(mat))
  expect_setequal(colnames(mat2), colnames(mat))

  # coerce to dgTMatrix so we can compare directly
  mat1 <- as(mat, "dgTMatrix")
  rlabs <- rownames(mat2)
  clabs <- colnames(mat2)
  expect_equal(mat1[rlabs, clabs], mat2[rlabs, clabs])
})

test_that("array dimensions can be transposed", {
  uri <- withr::local_tempdir("assay-matrix-transpose")
  mat <- create_sparse_matrix_with_string_dims(nrows = 10, ncols = 5)

  assaymat <- AssayMatrix$new(uri = uri, verbose = FALSE)
  assaymat$from_matrix(mat, index_cols = c("i", "j"), transpose = TRUE)
  expect_equal(assaymat$dimnames(), c("j", "i"))

  df2 <- assaymat$to_dataframe()
  # dim values match dim names
  expect_equal(unique(df2$j), colnames(mat))
  expect_setequal(unique(df2$i), rownames(mat))

  # recreated matrix is transposed
  mat2 <- assaymat$to_matrix()
  expect_equal(rownames(mat2), colnames(mat))
  expect_setequal(colnames(mat2), rownames(mat))
})

test_that("Incomplete queries can be completed via batching", {
  uri <- withr::local_tempdir("assay-matrix-batched")
  with_allocation_size_preference(5e5)

  nr <- 1e3
  nc <- 1e2
  set.seed(1)
  smat <- Matrix::rsparsematrix(
    nrow = nr,
    ncol = nc,
    density = 0.8,
    rand.x = function(n) as.integer(runif(n, min = 1, max = 100)),
    repr = "T"
  )
  dimnames(smat) <- list(paste0("i", seq_len(nr)), paste0("j", seq_len(nc)))

  assaymat <- AssayMatrix$new(uri = uri, verbose = TRUE)
  assaymat$from_matrix(smat, index_cols = c("i", "j"), value_col = "counts")

  df1 <- assaymat$to_dataframe(batch_mode = FALSE)
  df2 <- assaymat$to_dataframe(batch_mode = TRUE)
  expect_equal(dim(df1), dim(df2))
})

test_that("user allocation has returned to its default size", {
  expect_equal(tiledb::get_allocation_size_preference(), 10485760)
})
