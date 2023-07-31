test_that("matrixZeroBasedView", {
  sm <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = c(41, 42, 43))
  dm <- as.matrix(sm)
  for (mat1 in list(sm, dm)) {
    mat <- matrixZeroBasedView$new(mat1)

    # Test row and column indexing
    expect_equal(mat$take(0, 0)$get_one_based_matrix()[,,drop=T], 41)
    expect_equal(mat$take(1, 1)$get_one_based_matrix()[,,drop=T], 42)
    expect_equal(mat$take(2, 2)$get_one_based_matrix()[,,drop=T], 43)
    expect_equal(mat$take(0, 1)$get_one_based_matrix()[,,drop=T], 0)
    expect_equal(mat$take(2, 0)$get_one_based_matrix()[,,drop=T], 0)
    
    # Test row indexing only
    expect_equal(mat$take(i=1)$get_one_based_matrix()[,,drop=T], c(0, 42, 0))

    # Test column indexing only
    expect_equal(mat$take(j=1)$get_one_based_matrix()[,,drop=T], c(0, 42, 0))

    # Test vector slicing.
    slice <- mat$take(c(0, 2), c(0, 2))$get_one_based_matrix()
    expect_equal(dim(slice), c(2, 2))
    # slice itself is one-based!
    expect_equal(slice[1, ], c(41, 0))
    expect_equal(slice[, 1], c(41, 0))
    expect_equal(slice[2, ], c(0, 43))
    expect_equal(slice[, 2], c(0, 43))

    slice <- mat$take(i=c(0, 2))$get_one_based_matrix()[,,drop=TRUE]
    expect_equal(dim(slice), c(2, 3))
    expect_equal(slice[1, ], c(41, 0, 0))
    expect_equal(slice[2, ], c(0, 0, 43))

    slice <- mat$take(j=c(0, 2))$get_one_based_matrix()[,,drop=TRUE]
    expect_equal(dim(slice), c(3, 2))
    expect_equal(slice[, 1], c(41, 0, 0))
    expect_equal(slice[, 2], c(0, 0, 43))

    slice <- mat$take()$get_one_based_matrix()
    expect_equal(dim(slice), c(3, 3))
    expect_equal(slice[1, 1], 41)

    # Test misc properties
    expect_equal(mat$dim(), c(3, 3))
    expect_equal(mat$nrow(), 3)
    expect_equal(mat$ncol(), 3)

    # Test addition
    expect_equal(mat$sum(mat)$get_one_based_matrix(), mat1 + mat1)
    expect_error(mat$sum(1))
  }
})
