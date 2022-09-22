
test_that("conversion of dgTMatrix to COO data frame", {
  # create an ordered dgTMatrix by starting with an ordered dgCMatrix and
  # coercing that to a dgTMatrix
  omat <- as(create_sparse_matrix_with_string_dims(repr = "C"), "dgTMatrix")

  # create an unordered dgTMatrix
  umat <- create_sparse_matrix_with_string_dims(repr = "T")

  # verify coorindates are in column-major order
  expect_true(
    identical(
      omat@x,
      Filter(function(x) x != 0, as.numeric(omat))
    )
  )

  # verify coorindates are *not* in column-major order
  expect_false(
    identical(
      umat@x,
      Filter(function(x) x != 0, as.numeric(umat))
    )
  )

  # but materialized coordinates are identical
  expect_true(all(umat == omat))

  # verify round-tripping of omat to coo and back
  odf <- matrix_to_coo(omat)
  expect_true(is.data.frame(odf))

  ilabs <- unique(odf$i)
  jlabs <- unique(odf$j)
  expect_setequal(ilabs, rownames(omat))
  expect_setequal(jlabs, colnames(omat))

  omat2 <- dataframe_to_dgtmatrix(odf)[[1]]
  expect_identical(
    omat2[ilabs, jlabs],
    omat[ilabs, jlabs]
  )

  # verify round-tripping of umat to coo and back
  udf <- matrix_to_coo(umat)
  umat2 <- dataframe_to_dgtmatrix(udf)[[1]]
  expect_identical(
    umat2[ilabs, jlabs],
    umat[ilabs, jlabs]
  )

  # convert list of dgtMatrix objects with heterogenous coordinate ordering
  df <- matrix_to_coo(list(ordered = omat, unordered = umat))
  expect_identical(df$ordered, df$unordered)
})

test_that("conversion of a list dgTMatrix's to COO data frame", {
  mats <- list(
    SeuratObject::GetAssayData(pbmc_small, "counts"),
    SeuratObject::GetAssayData(pbmc_small, "data")
  )
  mats <- lapply(mats, FUN = as, Class = "dgTMatrix")

  df <- matrix_to_coo(mats)
  testthat::expect_true(is.data.frame(df))
  testthat::expect_equal(ncol(df), 4)

  ilabs <- unique(df$i)
  jlabs <- unique(df$j)

  mats2 <- dataframe_to_dgtmatrix(df, index_cols = c("i", "j"))
  expect_identical(
    mats[[1]][ilabs, jlabs],
    mats2[[1]][ilabs, jlabs]
  )
  expect_identical(
    mats[[2]][ilabs, jlabs],
    mats2[[2]][ilabs, jlabs]
  )

  mats[[3]] <- SeuratObject::GetAssayData(pbmc_small, "scale.data")

  mats[[3]] <- as(mats[[3]], "dgTMatrix")
  expect_error(
    matrix_to_coo(mats),
    "Matrix 1 and 3 are not layerable"
  )
})
