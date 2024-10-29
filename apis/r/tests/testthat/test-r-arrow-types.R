test_that("Arrow to R types: data type", {
  skip_if(!extended_tests())
  skip_if_not_installed("arrow")

  ints <- apply(
    expand.grid(c("", "u"), "int", c("8", "16", "32")),
    MARGIN = 1L,
    FUN = paste,
    collapse = ""
  )
  for (i in c(ints, "dictionary")) {
    f <- get(i, envir = asNamespace("arrow"))
    expect_type(rt <- r_type_from_arrow_type(f()), "character")
    expect_length(rt, 1L)
    expect_null(names(rt))
    expect_identical(
      rt,
      "integer",
      label = sprintf("r_type_from_arrow_type(arrow::%s())", i)
    )
  }

  dbls <- c("int64", "uint64", "date32", "timestamp", "float", "float32")
  for (i in dbls) {
    f <- get(i, envir = asNamespace("arrow"))
    expect_type(rt <- r_type_from_arrow_type(f()), "character")
    expect_length(rt, 1L)
    expect_null(names(rt))
    expect_identical(
      rt,
      "double",
      label = sprintf("r_type_from_arrow_type(arrow::%s())", i)
    )
  }

  for (i in c("bool", "boolean")) {
    f <- get(i, envir = asNamespace("arrow"))
    expect_type(rt <- r_type_from_arrow_type(f()), "character")
    expect_length(rt, 1L)
    expect_null(names(rt))
    expect_identical(
      rt,
      "logical",
      label = sprintf("r_type_from_arrow_type(arrow::%s())", i)
    )
  }
  for (i in c("utf8", "string", "large_utf8")) {
    f <- get(i, envir = asNamespace("arrow"))
    expect_type(rt <- r_type_from_arrow_type(f()), "character")
    expect_length(rt, 1L)
    expect_null(names(rt))
    expect_identical(
      rt,
      "character",
      label = sprintf("r_type_from_arrow_type(arrow::%s())", i)
    )
  }
})

test_that("Arrow to R types: field", {
  skip_if(!extended_tests())
  skip_if_not_installed("arrow")

  field <- arrow::field(name = random_name(), type = arrow::int8())
  expect_type(rt <- r_type_from_arrow_type(field), "character")
  expect_length(rt, 1L)
  expect_named(rt, field$name)
  expect_equivalent(rt, "integer")
})

test_that("Arrow to R types: schema", {
  skip_if(!extended_tests())

  asch <- create_arrow_schema()
  expect_type(rt <- r_type_from_arrow_type(asch), "character")
  expect_length(rt, length(asch))
  expect_named(rt, asch$names)
  for (fn in names(rt)) {
    et <- switch(
      EXPR = fn,
      int_column = "integer",
      soma_joinid = "double",
      float_column = "double",
      string_column = "character"
    )
    expect_equivalent(
      rt[fn],
      et,
      label = sprintf("r_type_from_arrow_type(schema[[%s]])", fn),
      expected.label = dQuote(et, FALSE)
    )
  }
})
