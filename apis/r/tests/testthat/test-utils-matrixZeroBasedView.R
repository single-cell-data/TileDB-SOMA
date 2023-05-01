test_that("matrixZeroBasedView", {
  sm <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = c(41, 42, 43))
  dm <- as.matrix(sm)
  for (mat1 in list(sm, dm)) {
    mat <- matrixZeroBasedView(mat1)

    # Test row and column indexing
    expect_equal(mat[0, 0], 41)
    expect_equal(mat[1, 1], 42)
    expect_equal(mat[2, 2], 43)
    expect_equal(mat[0, 1], 0)
    expect_equal(mat[2, 0], 0)

    # Test row indexing only
    expect_equal(mat[1, ], c(0, 42, 0))

    # Test column indexing only
    expect_equal(mat[, 1], c(0, 42, 0))

    # Test vector slicing.
    slice <- mat[c(0, 2), c(0, 2)]
    expect_equal(dim(slice), c(2, 2))
    # slice itself is one-based!
    expect_equal(slice[1, ], c(41, 0))
    expect_equal(slice[, 1], c(41, 0))
    expect_equal(slice[2, ], c(0, 43))
    expect_equal(slice[, 2], c(0, 43))

    slice <- mat[c(0, 2), ]
    expect_equal(dim(slice), c(2, 3))
    expect_equal(slice[1, ], c(41, 0, 0))
    expect_equal(slice[2, ], c(0, 0, 43))

    slice <- mat[, c(0, 2)]
    expect_equal(dim(slice), c(3, 2))
    expect_equal(slice[, 1], c(41, 0, 0))
    expect_equal(slice[, 2], c(0, 0, 43))

    slice <- mat[, ]
    expect_equal(dim(slice), c(3, 3))
    expect_equal(slice[1, 1], 41)

    # Test misc properties
    expect_equal(dim(mat), c(3, 3))
    expect_equal(nrow(mat), 3)
    expect_equal(ncol(mat), 3)

    # reject mutation
    rdo <- "matrixZeroBasedView is read-only"
    expect_error(mat[1, 1] <- 99, rdo)
    expect_error(mat[1, 1] <- NA, rdo)
    expect_error(mat[1, ] <- c(0, 99, 0), rdo)
    expect_error(mat[, 1] <- c(0, 99, 0), rdo)

    # reject arithmetic
    expect_error(mat + 1)
    expect_error(mat + mat)
  }
})
