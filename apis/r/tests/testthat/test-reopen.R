test_that("Test reopen works on arrays", {
  shape <- c(500L, 100L)
  for (cls in c("SOMADataFrame", "SOMASparseNDArray", "SOMADenseNDArray")) {
    uri <- withr::local_tempdir(paste("soma", cls, "array", "reopen", sep = "-"))
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
