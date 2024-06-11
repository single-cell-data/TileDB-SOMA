test_that("`reopen()` works on arrays", {
  shape <- c(500L, 100L)
  for (cls in c("SOMADataFrame", "SOMASparseNDArray", "SOMADenseNDArray")) {
    uri <- tempfile(pattern=paste("soma", cls, "reopen", sep = "-"))
    arr <- switch(
      EXPR = cls,
      SOMADataFrame = SOMADataFrameCreate(
        uri,
        schema = arrow::infer_schema(data.frame(
          soma_joinid = bit64::integer64(),
          int = integer()
        ))
      ),
      SOMASparseNDArray = SOMASparseNDArrayCreate(
        uri,
        type = arrow::int32(),
        shape = shape
      ),
      SOMADenseNDArray = SOMADenseNDArrayCreate(
        uri,
        type = arrow::int32(),
        shape = shape
      )
    )
    expect_s3_class(arr, cls)
    expect_identical(
      arr$mode(),
      "WRITE",
      info = sprintf("%sCreate() returns object open for 'WRITE'", cls)
    )
    expect_true(
      arr$is_open(),
      info = sprintf("%s is open when mode is 'WRITE'", cls)
    )

    lab <- sprintf("%s$reopen()", cls)
    is_open <- sprintf("%s$reopen() returns an object object", cls)

    # Test implicit WRITE -> READ
    expect_invisible(arr$reopen(), label = lab)
    expect_identical(
      arr$mode(),
      "READ",
      info = sprintf("%s$reopen() when mode is 'WRITE' reopens as 'READ'", cls)
    )
    expect_true(arr$is_open(), info = is_open)

    # Test implicit READ -> WRITE
    expect_invisible(arr$reopen(), label = lab)
    expect_identical(
      arr$mode(),
      "WRITE",
      info = sprintf("%s$reopen() when mode is 'READ' reopens as 'WRITE'", cls)
    )
    expect_true(arr$is_open(), info = is_open)

    # Test reopening in the same mode
    orig <- arr$mode()
    expect_invisible(arr$reopen(orig), label = lab)
    expect_identical(
      arr$mode(),
      orig,
      info = sprintf("%s$reopen() with an identical mode returns the same mode", cls)
    )
    expect_true(arr$is_open(), info = is_open)

    # Test reopen from close
    expect_no_condition(arr$close())
    expect_identical(arr$mode(), "CLOSED")
    expect_false(arr$is_open())

    for (mode in c("READ", "WRITE")) {
      expect_invisible(
        arr$reopen(mode),
        label = sprintf("%s$reopen('%s') from closed", cls, mode)
      )
      expect_identical(
        arr$mode(),
        mode,
        info = sprintf("%s$reopen('%s') returns a mode of '%s'", cls, mode, mode)
      )
      expect_true(
        arr$is_open(),
        info = sprintf("%s$reopen('%s') from CLOSED returns an open object", cls, mode)
      )
      expect_no_condition(arr$close())
    }

    arr$close()
    expect_error(arr$reopen())

    # Test assertions
    expect_error(arr$reopen("tomato"))
    expect_error(arr$reopen(TRUE))
    expect_error(arr$reopen(1L))

    arr$close()
  }
})

test_that("`reopen()` works on collections", {
  uri <- withr::local_tempdir("soma-collection-reopen")
  expect_s3_class(col <- SOMACollectionCreate(uri), "SOMACollection")
  expect_true(dir.exists(col$uri))

  expect_identical(col$mode(), "WRITE")
  expect_true(col$is_open())

  # Test implicit WRITE -> READ
  expect_invisible(col$reopen())
  expect_identical(col$mode(), "READ")
  expect_true(col$is_open())

  # Test implicit READ -> WRITE
  expect_invisible(col$reopen())
  expect_identical(col$mode(), "WRITE")
  expect_true(col$is_open())

  # Test reopening in the same mode
  orig <- col$mode()
  expect_invisible(col$reopen(orig))
  expect_identical(col$mode(), orig)
  expect_true(col$is_open())

  # Test reopen from close
  expect_no_condition(col$close())
  expect_identical(col$mode(), "CLOSED")
  expect_false(col$is_open())

  for (mode in c("READ", "WRITE")) {
    expect_invisible(col$reopen(mode))
    expect_identical(col$mode(), mode)
    expect_true(col$is_open())
    expect_no_condition(col$close())
  }

  col$close()
  expect_error(col$reopen())
})

test_that("`reopen()` works on nested collections", {
  shape <- c(500L, 100L)
  uri <- withr::local_tempdir("soma-nested-collection-reopen")
  expect_s3_class(col <- SOMACollectionCreate(uri), "SOMACollection")
  expect_true(dir.exists(col$uri))

  expect_identical(col$mode(), "WRITE")
  expect_true(col$is_open())

  col$add_new_collection(
    SOMACollectionCreate(file.path(col$uri, "collection")),
    "collection"
  )
  col$add_new_dense_ndarray("dense", type = arrow::int32(), shape = shape)
  col$add_new_sparse_ndarray("sparse", type = arrow::int32(), shape = shape)
  col$add_new_dataframe(
    "dataframe",
    schema = arrow::infer_schema(data.frame(
      soma_joinid = bit64::integer64(),
      int = integer()
    )),
    index_column_names = "soma_joinid"
  )
  expect_length(col$names(), 4L)

  for (i in col$names()) {
    expect_identical(
      col$get(i)$mode(),
      "WRITE",
      info = paste("subobject", sQuote(i), "has a mode of WRITE")
    )
    expect_true(
      col$get(i)$is_open(),
      info = paste("subobject", sQuote(i), "is open")
    )
  }

  # Test implicit WRITE -> READ
  expect_invisible(col$reopen())
  expect_identical(col$mode(), "READ")
  expect_true(col$is_open())

  for (i in col$names()) {
    expect_identical(
      col$get(i)$mode(),
      "READ",
      info = paste("subobject", sQuote(i), "has a mode of READ")
    )
    expect_true(
      col$get(i)$is_open(),
      info = paste("subobject", sQuote(i), "is open")
    )
  }

  # Test implicit READ -> WRITE
  expect_invisible(col$reopen())
  expect_identical(col$mode(), "WRITE")
  expect_true(col$is_open())

  for (i in col$names()) {
    expect_identical(
      col$get(i)$mode(),
      "WRITE",
      info = paste("subobject", sQuote(i), "has a mode of WRITE")
    )
    expect_true(
      col$get(i)$is_open(),
      info = paste("subobject", sQuote(i), "is open")
    )
  }

  # Test reopening in the same mode
  orig <- col$mode()
  expect_invisible(col$reopen(orig))
  expect_identical(col$mode(), orig)
  expect_true(col$is_open())

  for (i in col$names()) {
    expect_identical(
      col$get(i)$mode(),
      orig,
      info = paste("subobject", sQuote(i), "has a mode of", orig)
    )
    expect_true(
      col$get(i)$is_open(),
      info = paste("subobject", sQuote(i), "is open")
    )
  }

  # Test reopen from close
  expect_no_condition(col$close())
  expect_identical(col$mode(), "CLOSED")
  expect_false(col$is_open())

  # Cannot test that `close()` affects subobjects as
  # the collection is closed

  # for (i in col$names()) {
  #   expect_identical(
  #     col$get(i)$mode(),
  #     "CLOSED",
  #     info = paste("subobject", sQuote(i), "has a mode of CLOSED")
  #   )
  #   expect_false(
  #     col$get(i)$is_open(),
  #     info = paste("subobject", sQuote(i), "is closed")
  #   )
  # }

  for (mode in c("READ", "WRITE")) {
    expect_invisible(col$reopen(mode))
    expect_identical(col$mode(), mode)
    expect_true(col$is_open())
    for (i in col$names()) {
      expect_identical(
        col$get(i)$mode(),
        mode,
        info = paste("subobject", sQuote(i), "has a mode of", mode)
      )
      expect_true(
        col$get(i)$is_open(),
        info = paste("subobject", sQuote(i), "is open")
      )
    }
    expect_no_condition(col$close())
  }

  col$close()
  expect_error(col$reopen())
})

test_that("`reopen()` works on SOMAMeasurements", {
  uri <- withr::local_tempdir("soma-measurement-reopen")
  expect_s3_class(ms <- SOMAMeasurementCreate(uri), "SOMAMeasurement")
  expect_true(dir.exists(ms$uri))

  expect_identical(ms$mode(), "WRITE")
  expect_true(ms$is_open())

  ms$var <- create_and_populate_var(file.path(ms$uri, "var"))
  ms$X <- SOMACollectionCreate(file.path(ms$uri, "X"))

  ms$X$set(
    create_and_populate_sparse_nd_array(file.path(ms$X$uri, "RNA")),
    name = "RNA"
  )

  # Test implicit WRITE -> READ
  expect_invisible(ms$reopen())
  expect_identical(ms$mode(), "READ")
  expect_true(ms$is_open())

  expect_identical(ms$var$mode(), "READ")
  expect_true(ms$var$is_open())

  expect_identical(ms$X$mode(), "READ")
  expect_true(ms$X$is_open())

  # Test implicit READ -> WRITE
  expect_invisible(ms$reopen())
  expect_identical(ms$mode(), "WRITE")
  expect_true(ms$is_open())

  expect_identical(ms$var$mode(), "WRITE")
  expect_true(ms$var$is_open())

  expect_identical(ms$X$mode(), "WRITE")
  expect_true(ms$X$is_open())

  # Test reopening in the same mode
  orig <- ms$mode()
  expect_invisible(ms$reopen(orig))
  expect_identical(ms$mode(), orig)
  expect_true(ms$is_open())

  expect_identical(ms$var$mode(), orig)
  expect_true(ms$var$is_open())

  expect_identical(ms$X$mode(), orig)
  expect_true(ms$X$is_open())

  # Test reopen from close
  expect_no_condition(ms$close())
  expect_identical(ms$mode(), "CLOSED")
  expect_false(ms$is_open())

  # Cannot test that `close()` affects subobjects as
  # the measurement is closed

  # expect_identical(ms$var$mode(), "CLOSED")
  # expect_false(ms$var$is_open())
  #
  # expect_identical(ms$X$mode(), "CLOSED")
  # expect_false(ms$X$is_open())

  for (mode in c("READ", "WRITE")) {
    expect_invisible(ms$reopen(mode))
    expect_identical(ms$mode(), mode)
    expect_true(ms$is_open())
    expect_identical(ms$var$mode(), mode)
    expect_true(ms$var$is_open())
    expect_identical(ms$X$mode(), mode)
    expect_true(ms$X$is_open())
    expect_no_condition(ms$close())
  }

  ms$close()
  expect_error(ms$reopen())
})

test_that("`reopen()` works on SOMAExperiments", {
  uri <- withr::local_tempdir("soma-experiment-reopen")
  exp <- create_and_populate_experiment(
    uri = uri,
    n_obs = 20L,
    n_var = 10L,
    X_layer_names = "counts",
    mode = "WRITE"
  )
  expect_identical(exp$mode(), "WRITE")
  expect_true(exp$is_open())

  # Test implicit WRITE -> READ
  expect_invisible(exp$reopen())
  expect_identical(exp$mode(), "READ")
  expect_true(exp$is_open())

  expect_identical(exp$obs$mode(), "READ")
  expect_true(exp$obs$is_open())

  expect_identical(exp$ms$mode(), "READ")
  expect_true(exp$ms$is_open())

  # Test implicit READ -> WRITE
  expect_invisible(exp$reopen())
  expect_identical(exp$mode(), "WRITE")
  expect_true(exp$is_open())

  expect_identical(exp$obs$mode(), "WRITE")
  expect_true(exp$obs$is_open())

  expect_identical(exp$ms$mode(), "WRITE")
  expect_true(exp$ms$is_open())

  # Test reopening in the same mode
  orig <- exp$mode()
  expect_invisible(exp$reopen(orig))
  expect_identical(exp$mode(), orig)
  expect_true(exp$is_open())

  expect_identical(exp$obs$mode(), orig)
  expect_true(exp$obs$is_open())

  expect_identical(exp$ms$mode(), orig)
  expect_true(exp$ms$is_open())

  # Test reopen from close
  expect_no_condition(exp$close())
  expect_identical(exp$mode(), "CLOSED")
  expect_false(exp$is_open())

  # Cannot test that `close()` affects subobjects as
  # the measurement is closed

  # expect_identical(exp$obs$mode(), "CLOSED")
  # expect_false(exp$obs$is_open())
  #
  # expect_identical(exp$ms$mode(), "CLOSED")
  # expect_false(exp$ms$is_open())

  for (mode in c("READ", "WRITE")) {
    expect_invisible(exp$reopen(mode))
    expect_identical(exp$mode(), mode)
    expect_true(exp$is_open())
    expect_identical(exp$obs$mode(), mode)
    expect_true(exp$obs$is_open())
    expect_identical(exp$ms$mode(), mode)
    expect_true(exp$ms$is_open())
    expect_no_condition(exp$close())
  }

  exp$close()
  expect_error(exp$reopen())

})
