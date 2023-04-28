# This sparseMatrix subclass overrides the [ accessor to take zero-based coordinates
# instead of the usual R one-based convention. This is needed for SOMASparseNDArray
# with soma_joinid dimensions, since soma_joinid can be zero (and it's not appropriate
# to implicitly increment it, since it's an id & not an offset).
# See discussion: https://github.com/single-cell-data/TileDB-SOMA/issues/1232

# We'll do some weird stuff to subclass -each- of these sparseMatrix implementations:
sparseMatrixClasses <- c(
  "dgCMatrix", "dgTMatrix", "dgRMatrix",
  "lgCMatrix", "lgTMatrix", "lgRMatrix"
)

# factory to declare the subclass for one of the above superclasses
sparseMatrixZeroBasedClassDeclaration <- function(superclass) {
  stopifnot(superclass %in% sparseMatrixClasses)
  subclass <- paste(superclass, "ZeroBased", sep = "")
  setClass(subclass, contains = getClass(superclass, where = "Matrix"))
  setMethod(
    "[", signature(
      x = subclass,
      i = "index", j = "index", drop = "logical"
    ),
    function(x, i, j, drop) {
      as(x, superclass)[i + 1, j + 1, drop = drop]
      # The returned vector/matrix will be one-based! While we could add a zero-based
      # wrapper to it here, that would open a can of worms about all possible matrix-
      # valued operations on x (e.g. x+1, x!=0, ...)
    }
  )
  setMethod(
    "[", signature(
      x = subclass,
      i = "index", j = "missing", drop = "logical"
    ),
    function(x, i, j, drop) {
      stopifnot(missing(j))
      as(x, superclass)[i + 1, , drop = drop]
    }
  )
  setMethod(
    "[", signature(
      x = subclass,
      i = "missing", j = "index", drop = "logical"
    ),
    function(x, i, j, drop) {
      stopifnot(missing(i))
      as(x, superclass)[, j + 1, drop = drop]
    }
  )
  # for consistency: x[,] has to return one-based too
  setMethod(
    "[", signature(x = subclass, i = "missing", j = "missing", drop = "logical"),
    function(x, i, j, drop) {
      stopifnot(missing(i) && missing(j))
      as(x, superclass)
    }
  )
  # disable mutation
  for (i_ty in c("index", "missing")) {
    for (j_ty in c("index", "missing")) {
      for (value_ty in c("replValue", "sparseVector")) {
        setReplaceMethod(
          "[", signature(x = subclass, i = i_ty, j = j_ty, value = value_ty),
          function(x, i, j, ..., value) {
            stop("sparseMatrixZeroBased is a read-only view.")
          }
        )
      }
    }
  }
}

# operate the factory
for (superclass in sparseMatrixClasses) {
  sparseMatrixZeroBasedClassDeclaration(superclass)
}

# Given a (one-based) sparseMatrix x, return the zero-based view of it.
sparseMatrixZeroBasedView <- function(x) {
  stopifnot(class(x) %in% sparseMatrixClasses)
  new(paste(class(x), "ZeroBased", sep = ""), x)
}

# as Matrix::sparseMatrix(..., index1 = FALSE), but return the zero-based view
sparseMatrixZeroBased <- function(...) {
  sparseMatrixZeroBasedView(sparseMatrix(..., index1 = FALSE))
}
