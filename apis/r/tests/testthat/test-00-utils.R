test_that("validate read coords", {
  # NULL is a valid value
  expect_equal(
    validate_read_coords(NULL),
    NULL
  )

  # bare vector is converted to list
  expect_equal(
    validate_read_coords(1:10),
    list(1:10)
  )

  expect_equal(
    validate_read_coords(list(1:10)),
    list(1:10)
  )

  # list of multiple coordinate vectors must be named
  expect_error(
    validate_read_coords(list(1:10, 1:10)),
    "'coords' must be a named list to query multiple dimensions"
  )

  expect_equal(
    validate_read_coords(list(int_column = 1:10, float_column = 1:10)),
    list(int_column = 1:10, float_column = 1:10)
  )
})

test_that("validate read coords with dimension names", {
  # assume vector or unnamed list of length 1 corresponds to first dimension
  expect_equal(
    validate_read_coords(1:10, dimnames = "int_column"),
    list(int_column = 1:10)
  )

  # list of named coordinates must match provided dimension names
  expect_error(
    validate_read_coords(list(int_column = 1:10, float_column = 1:10), c("int_column", "string_column")),
    "names of 'coords' must correspond to dimension names"
  )

  expect_error(
    validate_read_coords(list(int_column = 1:10, float_column = 1:10), c("int_column")),
    "names of 'coords' must correspond to dimension names"
  )
})

test_that("validate read coords with dimension names and schema", {
  asch <- create_arrow_schema()

  # if schema is provided, dimnames must be provided too
  expect_error(
    validate_read_coords(1:10, schema = asch),
    "'dimnames' must be provided with a 'schema'"
  )

  # assume vector or unnamed list of length 1 corresponds to first dimension
  expect_equal(
    validate_read_coords(1:10, dimnames = "int_column", schema = asch),
    list(int_column = 1:10)
  )

  # integer coordinates corresponding to int64 dimensions are cast to int64
  expect_equal(
    validate_read_coords(1:10, dimnames = "soma_joinid", schema = asch),
    list(soma_joinid = as.integer64(1:10))
  )

  # casting is selective and only applies to int64 dimensions
  test_coords <- validate_read_coords(
    coords = list(int_column = 1:10, soma_joinid = 1:10),
    dimnames = c("int_column", "soma_joinid"),
    schema = asch
  )

  expect_equal(test_coords$int_column, 1:10)
  expect_equal(test_coords$soma_joinid, as.integer64(1:10))
})

test_that("half-named lists are not treated as named", {
  expect_true(is_named_list(list(a = 1, b = 2)))
  expect_false(is_named_list(list(a = 1, 2)))
  expect_false(is_named_list(list(1, 2)))
})

test_that("is_integerish: default", {
  expect_true(.is_integerish(vector("integer")))
  expect_true(.is_integerish(vector("numeric")))
  types <- c("logical", "complex", "character", "expression", "list", "raw")
  for (tt in types) {
    expect_false(
      .is_integerish(vector(tt)),
      label = sprintf(".is_integerish(vector('%s'))", tt)
    )
  }
})

test_that("is_integerish: integer64", {
  # Basic tests
  for (n in 0:3) {
    expect_true(
      .is_integerish(bit64::integer64(length = n)),
      label = sprintf(".is_integerish(integer64(length = %i))", n)
    )
    expect_true(
      .is_integerish(bit64::integer64(length = n), n = n),
      label = sprintf(".is_integerish(integer64(length = %i), n = %i)", n, n)
    )
    expect_false(
      .is_integerish(bit64::integer64(length = n), n = n + 1L),
      label = sprintf(".is_integerish(integer64(length = %i), n = %i)", n, n + 1L)
    )
  }

  # Test finiteness
  expect_true(.is_integerish(bit64::NA_integer64_))
  expect_true(.is_integerish(bit64::NA_integer64_, finite = FALSE))
  expect_false(.is_integerish(bit64::integer64(), finite = FALSE))
  expect_false(.is_integerish(bit64::NA_integer64_, finite = TRUE))

  # Test large number
  expect_true(.is_integerish(bit64::as.integer64((2^31) + 1L)))
})

test_that("is_integerish: arrow::DataType", {
  ints <- paste0("int", c(8, 16, 32, 64))
  for (it in c(ints, paste0("u", ints))) {
    f <- get(it, envir = asNamespace("arrow"))
    expect_true(
      .is_integerish(f()),
      label = sprintf(".is_integerish(arrow::%s())", it)
    )
  }

  types <- c(
    paste0("float", c("", 16, 32, 64)),
    paste0("bool", c("", "ean")),
    "string",
    paste0(c("", "large_"), "utf8"),
    paste0("date", c(32, 64)),
    paste0("time", c(32, 64, "stamp"))
  )
  for (at in types) {
    f <- get(at, envir = asNamespace("arrow"))
    expect_false(
      .is_integerish(f()),
      label = sprintf(".is_integerish(arrow::%s())", at)
    )
  }
})

test_that("is_integerish: arrow::Field", {
  sch <- create_arrow_schema()
  for (i in names(sch)) {
    label <- sprintf(".is_integerish(sch[['%s']])", i)
    switch(
      EXPR = i,
      int_column = ,
      soma_joinid = expect_true(.is_integerish(sch[[i]]), label = label),
      expect_false(.is_integerish(sch[[i]]), label = label)
    )
  }
})

test_that("is_integerish: arrow::Arrays", {
  tbl <- create_arrow_table()
  for (i in names(tbl)) {
    label <- sprintf(".is_integerish(tbl[['%s']])", i)
    switch(
      EXPR = i,
      int_column = ,
      soma_joinid = expect_true(.is_integerish(tbl[[i]]), label = label),
      expect_false(.is_integerish(tbl[[i]]), label = label)
    )
  }
})
