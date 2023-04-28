test_that("sparseMatrixZeroBased", {
  for (repr in c("C", "T", "R")) {
    mat <- sparseMatrixZeroBased(
      i = c(0, 1, 2), j = c(0, 1, 2), x = c(41, 42, 43),
      repr = repr
    )

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

    slice <- mat[,]
    expect_equal(dim(slice), c(3,3))
    expect_equal(slice[1,1], 41)

    # test rowSums, colSums, which
    expect_equal(Matrix::rowSums(mat), c(41, 42, 43))
    expect_equal(Matrix::colSums(mat), c(41, 42, 43))
    expect_equal(Matrix::which(mat != 0), c(1, 5, 9)) # one-based!

    # reject mutation attempts
    rdonly <- "sparseMatrixZeroBased is a read-only view."
    expect_error(mat[1, 1] <- 99, rdonly)
    expect_error(mat[1, 1] <- NA, rdonly)
    expect_error(mat[1, ] <- c(0, 99, 0), rdonly)
    expect_error(mat[, 1] <- c(0, 99, 0), rdonly)
  }
})
