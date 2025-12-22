test_that("SOMADataFrame shape", {
  asch <- create_arrow_schema()

  index_column_name_choices <- list(
    "soma_joinid",
    c("soma_joinid", "int_column"),
    c("soma_joinid", "string_column"),
    c("string_column", "int_column"), # duplicate intentional to match `domain_at_create_choices`
    c("string_column", "int_column")
  )

  domain_at_create_choices <- list(
    list(soma_joinid = c(0, 999)),
    list(soma_joinid = c(0, 999), int_column = c(-10000, 10000)),
    list(soma_joinid = c(0, 999), string_column = NULL),
    list(string_column = NULL, int_column = c(-10000, 10000)),
    list(string_column = c("", ""), int_column = c(-10000, 10000))
  )

  # Check the test configs themselves to make sure someone (ahem, me)
  # didn't edit one without forgetting to edit the other
  expect_equal(
    length(index_column_name_choices),
    length(domain_at_create_choices)
  )

  for (i in seq_along(index_column_name_choices)) {
    index_column_names <- index_column_name_choices[[i]]

    uri <- withr::local_tempdir("soma-dataframe-shape")

    # Create
    if (dir.exists(uri)) {
      unlink(uri, recursive = TRUE)
    }

    domain_for_create <- domain_at_create_choices[[i]]

    sdf <- SOMADataFrameCreate(
      uri,
      asch,
      index_column_names = index_column_names,
      domain = domain_for_create
    )

    expect_true(sdf$exists())
    expect_true(dir.exists(uri))

    # Write
    tbl0 <- arrow::arrow_table(
      int_column = 1L:4L,
      soma_joinid = 1L:4L,
      float_column = 1.1:4.1,
      string_column = c("apple", "ball", "cat", "dog"),
      schema = asch
    )

    sdf$write(tbl0)
    sdf$close()

    sdf <- SOMADataFrameOpen(uri)

    # Check shape and maxshape et al.
    expect_true(sdf$tiledbsoma_has_upgraded_domain())
    expect_error(sdf$shape(), class = "notYetImplementedError")
    expect_error(sdf$maxshape(), class = "notYetImplementedError")

    # Not implemented this way per
    # https://github.com/single-cell-data/TileDB-SOMA/pull/2953#discussion_r1746125089
    # sjid_shape <- sdf$.maybe_soma_joinid_shape()
    # sjid_maxshape <- sdf$.maybe_soma_joinid_maxshape()
    soma_context_handle <- create_soma_context()
    sjid_shape <- maybe_soma_joinid_shape(sdf$uri, soma_context_handle)
    sjid_maxshape <- maybe_soma_joinid_maxshape(sdf$uri, soma_context_handle)

    if ("soma_joinid" %in% index_column_names) {
      # More testing to come on
      # https://github.com/single-cell-data/TileDB-SOMA/issues/2407
      expect_false(rlang::is_na(sjid_shape))
      expect_false(rlang::is_na(sjid_maxshape))
    } else {
      expect_true(rlang::is_na(sjid_shape))
      expect_true(rlang::is_na(sjid_maxshape))
    }

    # Check has_upgraded_domain
    expect_true(sdf$tiledbsoma_has_upgraded_domain())

    # Check domain and maxdomain
    dom <- sdf$domain()
    mxd <- sdf$maxdomain()

    # First check names
    expect_equal(names(dom), index_column_names)
    expect_equal(names(mxd), index_column_names)

    # Then check all slots are pairs
    for (name in names(dom)) {
      expect_length(dom[[name]], 2L)
      expect_length(mxd[[name]], 2L)
    }

    # Then check contents

    # Old shape/domainishes (without core current domain) for non-string dims:
    # * There is no core current domain
    # * Expect domain == maxdomain
    # * If they asked for NULL: both should be huge (near min/max for datatype)
    # * If they asked for something specific: they should get it
    #
    # New shape/domainishes (with core current domain) for non-string dims:
    # * Maxdomain should be huge (near min/max for datatype)
    # * If they asked for NULL: domain should be the same as maxdomain
    # * If they asked for a specific domain: they should get it
    #
    # Old shape/domainishes (without core current domain) for string dims:
    # * There is no core current domain
    # * Expect domain == maxdomain
    # * Core domain for strings is always ("", "")
    #
    # New shape/domainishes (with core current domain) for string dims:
    # * Core domain (soma maxdomain) for strings is always ("", "")
    # * Core current domain (soma domain) for strings:
    #   o If they asked for NULL: expect ("", "")
    #   o If they asked for something specific: they should get it

    if ("soma_joinid" %in% index_column_names) {
      sjid_dom <- dom[["soma_joinid"]]
      sjid_mxd <- mxd[["soma_joinid"]]
      sjid_dfc <- domain_for_create[["soma_joinid"]]

      # Not: expect_equal(sjid_dom, bit64::as.integer64(sjid_dfc)) The
      # soma_joinid dim is always of type int64.  Everything coming back
      # from libtiledbsoma, through C nanoarrow, through the R arrow
      # package, to Arrow RecordBatch, holds true to that. But the final
      # as.list() converts the domain to regular integer. This is a feature
      # TBH: suppressable with `op <- options(arrow.int64_downcast =
      # FALSE)`. The maxdomainis likely to be in the 2**63 range
      # but the domain is likely to be ordinary-sized numbers in the
      # thousands or millions. Users are likely to prefer these
      # being downcast to regular R integers.
      expect_equal(sjid_dom, sjid_dfc)
    }

    if ("int_column" %in% index_column_names) {
      int_dom <- dom[["int_column"]]
      int_mxd <- mxd[["int_column"]]
      int_dfc <- domain_for_create[["int_column"]]

      expect_equal(int_dom, int_dfc)

      expect_true(int_mxd[[1]] < -2000000000)
      expect_true(int_mxd[[2]] > 2000000000)
    }

    if ("string_column" %in% index_column_names) {
      str_dom <- dom[["string_column"]]
      str_mxd <- mxd[["string_column"]]
      str_dfc <- domain_for_create[["string_column"]]

      if (is.null(str_dfc)) {
        expect_equal(str_dom, c("", ""))
      } else {
        expect_equal(str_dom, str_dfc)
      }
      expect_equal(str_mxd, c("", ""))
    }

    sdf$close()

    # Test resize for dataframes (more general upgrade_domain to be tested
    # separately -- see https://github.com/single-cell-data/TileDB-SOMA/issues/2407)
    has_soma_joinid_dim <- "soma_joinid" %in% index_column_names
    sjid_dfc <- domain_for_create[["soma_joinid"]]

    # Test resize down
    new_shape <- 0
    sdf <- SOMADataFrameOpen(uri, "WRITE")
    if (has_soma_joinid_dim) {
      # It's an error to downsize
      expect_error(sdf$tiledbsoma_resize_soma_joinid_shape(new_shape))
    } else {
      # There is no problem when soma_joinid is not a dim --
      # sdf$tiledbsoma_resize_soma_joinid_shape is a no-op in that case
      expect_no_condition(sdf$tiledbsoma_resize_soma_joinid_shape(new_shape))
    }
    sdf$close()

    # Make sure the failed resize really didn't change the shape
    if (has_soma_joinid_dim) {
      sdf <- SOMADataFrameOpen(uri, "READ")
      expect_equal(sdf$domain()[["soma_joinid"]], sjid_dfc)
      sdf$close()
    }

    # Test writes out of bounds, before resize
    old_shape <- 100
    if (has_soma_joinid_dim) {
      old_shape <- domain_for_create[["soma_joinid"]][[2]] + 1 + 100
    }
    new_shape <- old_shape + 100

    tbl1 <- arrow::arrow_table(
      int_column = 5L:8L,
      soma_joinid = (old_shape + 1L):(old_shape + 4L),
      float_column = 5.1:8.1,
      string_column = c("egg", "flag", "geese", "hay"),
      schema = asch
    )

    sdf <- SOMADataFrameOpen(uri, "WRITE")
    if (has_soma_joinid_dim) {
      expect_error(sdf$write(tbl1))
    } else {
      expect_no_condition(sdf$write(tbl1))
    }
    sdf$close()

    # Test resize
    sdf <- SOMADataFrameOpen(uri, "WRITE")
    sdf$tiledbsoma_resize_soma_joinid_shape(new_shape)
    sdf$close()

    # Test writes out of old bounds, within new bounds, after resize
    sdf <- SOMADataFrameOpen(uri, "WRITE")
    expect_no_condition(sdf$write(tbl1))
    sdf$close()

    # To do: test readback

    rm(tbl1)

    rm(sdf, tbl0)

    gc()
  }

  # Test `domain` assertions
  uri <- tempfile()

  # `domain` must be `NULL` or a list
  expect_error(SOMADataFrameCreate(
    uri,
    schema = asch,
    index_column_names = "soma_joinid",
    domain = NA
  ))
  expect_error(SOMADataFrameCreate(
    uri,
    schema = asch,
    index_column_names = "soma_joinid",
    domain = 1L
  ))
  expect_error(SOMADataFrameCreate(
    uri,
    schema = asch,
    index_column_names = "soma_joinid",
    domain = bit64::as.integer64(c(0L, 99L))
  ))
  # `domain` may not be an empty list
  expect_error(SOMADataFrameCreate(
    uri,
    schema = asch,
    index_column_names = "soma_joinid",
    domain = list()
  ))
  # `domain` must be named
  expect_error(SOMADataFrameCreate(
    uri,
    schema = asch,
    index_column_names = "soma_joinid",
    domain = list(bit64::as.integer64(c(0L, 99L)))
  ))
  # `domain` must be a list of two-length atomics
  expect_error(SOMADataFrameCreate(
    uri,
    schema = asch,
    index_column_names = "soma_joinid",
    domain = list(
      soma_joinid = list(bit64::as.integer64(0L), bit64::as.integer64(99L))
    )
  ))
  expect_error(SOMADataFrameCreate(
    uri,
    schema = asch,
    index_column_names = "soma_joinid",
    domain = list(soma_joinid = NA)
  ))
  expect_error(SOMADataFrameCreate(
    uri,
    schema = asch,
    index_column_names = "soma_joinid",
    domain = list(soma_joinid = numeric())
  ))
  expect_error(SOMADataFrameCreate(
    uri,
    schema = asch,
    index_column_names = "soma_joinid",
    domain = list(soma_joinid = numeric(length = 3L))
  ))
  expect_error(SOMADataFrameCreate(
    uri,
    schema = asch,
    index_column_names = c("soma_joinid", "int_column"),
    domain = list(soma_joinid = NULL, int_column = data.frame())
  ))
  # `names(domain)` must be identical to `index_column_names`
  expect_error(SOMADataFrameCreate(
    uri,
    schema = asch,
    index_column_names = "soma_joinid",
    domain = list(soma_joinid = NULL, int_column = NULL)
  ))
  expect_error(SOMADataFrameCreate(
    uri,
    schema = asch,
    index_column_names = c("soma_joinid", "int_column"),
    domain = list(soma_joinid = NULL)
  ))
})

test_that("SOMADataFrame domain mods", {
  uri <- withr::local_tempdir("soma-dataframe-domain-mods")

  schema = arrow::schema(
    arrow::field("soma_joinid", arrow::int64()),
    arrow::field("mystring", arrow::string()),
    arrow::field("myint", arrow::int16()),
    arrow::field("myfloat", arrow::float32()),
    arrow::field("mybool", arrow::bool()) # not supported as an index type
  )

  index_column_names <- c("soma_joinid", "mystring", "myint", "myfloat")

  domain_for_create <- list(
    soma_joinid = c(0, 3),
    mystring = NULL,
    myint = c(20, 50),
    myfloat = c(0.0, 6.0)
  )

  table <- arrow::arrow_table(
    soma_joinid = 0L:3L,
    mystring = c("a", "b", "a", "b"),
    myint = c(20, 30, 40, 50),
    myfloat = c(1.0, 2.5, 4.0, 5.5),
    mybool = c(TRUE, FALSE, TRUE, TRUE),
    schema = schema
  )

  sdf <- SOMADataFrameCreate(
    uri,
    schema = schema,
    index_column_names = index_column_names,
    domain = domain_for_create
  )
  sdf$write(table)
  sdf$close()

  sdf <- SOMADataFrameOpen(uri, "WRITE")

  # Check "expand" to same
  new_domain <- list(
    soma_joinid = c(0, 3),
    mystring = NULL,
    myint = c(20, 50),
    myfloat = c(0.0, 6.0)
  )

  # -- first check dry run
  expect_no_condition(sdf$change_domain(domain_for_create, check_only = TRUE))
  sdf$close()

  check <- list(
    soma_joinid = c(0, 3),
    mystring = c("", ""), # this is how it reads back
    myint = c(20, 50),
    myfloat = c(0.0, 6.0)
  )
  expect_equal(sdf$domain(), check)
  sdf$close()

  sdf <- SOMADataFrameOpen(uri, "WRITE")
  expect_no_condition(sdf$change_domain(new_domain))

  # Shrink
  new_domain <- list(
    soma_joinid = c(0, 2),
    mystring = NULL,
    myint = c(20, 50),
    myfloat = c(0.0, 6.0)
  )
  expect_error(sdf$change_domain(new_domain))

  new_domain <- list(
    soma_joinid = c(0, 3),
    mystring = NULL,
    myint = c(20, 40),
    myfloat = c(0.0, 6.0)
  )
  expect_error(sdf$change_domain(new_domain))

  new_domain <- list(
    soma_joinid = c(0, 3),
    mystring = NULL,
    myint = c(20, 50),
    myfloat = c(2.0, 6.0)
  )
  expect_error(sdf$change_domain(new_domain))

  # String domain cannot be specified
  new_domain <- list(
    soma_joinid = c(0, 3),
    mystring = c("a", "z"),
    myint = c(20, 50),
    myfloat = c(0.0, 6.0)
  )
  expect_error(sdf$change_domain(new_domain))

  # All clear

  # String domain cannot be specified
  new_domain <- list(
    soma_joinid = c(0, 8),
    mystring = c("", ""),
    myint = c(20, 50),
    myfloat = c(0.0, 10.0)
  )
  expect_no_condition(sdf$change_domain(new_domain))
  sdf$close()

  # Check for success
  sdf <- SOMADataFrameOpen(uri, "READ")
  dom <- sdf$domain()
  expect_equal(sdf$domain(), new_domain)
  sdf$close()
})

test_that("SOMASparseNDArray shape", {
  uri <- withr::local_tempdir("soma-sparse-ndarray-shape")
  asch <- create_arrow_schema()

  element_type_choices <- list(arrow::float32(), arrow::int16())
  arg_shape <- c(100, 200)
  for (element_type in element_type_choices) {
    if (dir.exists(uri)) {
      unlink(uri, recursive = TRUE)
    }
    ndarray <- SOMASparseNDArrayCreate(uri, element_type, shape = arg_shape)
    ndarray$close()

    ndarray <- SOMASparseNDArrayOpen(uri)
    expect_equal(ndarray$ndim(), length(arg_shape))

    expect_equal(ndarray$shape(), bit64::as.integer64(arg_shape))

    # More generally after current-domain support:
    readback_shape <- ndarray$shape()
    readback_maxshape <- ndarray$maxshape()
    expect_equal(length(readback_shape), length(readback_maxshape))
    expect_true(all(readback_shape < readback_maxshape))

    ndarray$close()

    # Test write in bounds
    ndarray <- SOMASparseNDArrayOpen(uri, "WRITE")
    soma_dim_0 <- c(2, 3)
    soma_dim_1 <- c(4, 5)
    soma_data <- c(60, 70)
    sm <- Matrix::sparseMatrix(i = soma_dim_0, j = soma_dim_1, x = soma_data)
    ndarray$write(sm)
    ndarray$close()

    ndarray <- SOMASparseNDArrayOpen(uri)
    ned <- ndarray$non_empty_domain(max_only = TRUE)
    # expect_equal(ned, c(2,4))
    expect_equal(as.integer(ned), as.integer(c(2, 4)))

    # Test reads out of bounds
    coords <- list(bit64::as.integer64(c(1, 2)), bit64::as.integer64(c(3, 4)))
    expect_no_error(x <- ndarray$read(coords = coords)$tables()$concat())

    coords <- list(
      bit64::as.integer64(c(101, 202)),
      bit64::as.integer64(c(3, 4))
    )
    expect_error(x <- ndarray$read(coords = coords)$tables()$concat())

    ndarray$close()

    ndarray <- SOMASparseNDArrayOpen(uri, "WRITE")

    # Test resize down
    new_shape <- c(50, 60)
    expect_error(ndarray$resize(new_shape))

    # Test writes out of old bounds
    soma_dim_0 <- c(200, 300)
    soma_dim_1 <- c(400, 500)
    soma_data <- c(6000, 7000)
    sm <- Matrix::sparseMatrix(i = soma_dim_0, j = soma_dim_1, x = soma_data)
    expect_error(ndarray$write(sm))

    # Test resize up
    new_shape <- c(500, 600)
    #### expect_no_error(ndarray$resize(new_shape))
    ndarray$resize(new_shape)

    # Test writes within new bounds
    soma_dim_0 <- c(200, 300)
    soma_dim_1 <- c(400, 500)
    soma_data <- c(6000, 7000)
    sm <- Matrix::sparseMatrix(i = soma_dim_0, j = soma_dim_1, x = soma_data)
    expect_no_error(ndarray$write(sm))
    ndarray$close()

    ndarray <- SOMASparseNDArrayOpen(uri)
    coords <- list(
      bit64::as.integer64(c(101, 202)),
      bit64::as.integer64(c(3, 4))
    )
    expect_no_error(x <- ndarray$read(coords = coords)$tables()$concat())
    ndarray$close()

    rm(ndarray)
    gc()
  }
})

test_that("SOMADenseNDArray shape", {
  skip_if(!extended_tests())
  skip_on_cran()
  uri <- withr::local_tempdir("soma-dense-ndarray-shape")
  asch <- create_arrow_schema()

  element_type_choices <- list(arrow::float32(), arrow::int16())
  arg_shape <- c(100, 200)
  for (element_type in element_type_choices) {
    if (dir.exists(uri)) {
      unlink(uri, recursive = TRUE)
    }
    ndarray <- SOMADenseNDArrayCreate(uri, element_type, shape = arg_shape)
    ndarray$close()

    ndarray <- SOMADenseNDArrayOpen(uri)

    expect_equal(ndarray$ndim(), length(arg_shape))

    expect_equal(ndarray$shape(), bit64::as.integer64(arg_shape))

    # More generally after current-domain support:
    readback_shape <- ndarray$shape()
    readback_maxshape <- ndarray$maxshape()
    expect_equal(length(readback_shape), length(readback_maxshape))

    expect_true(all(readback_shape < readback_maxshape))

    ndarray$close()

    # Test write in bounds
    ndarray <- SOMADenseNDArrayOpen(uri, "WRITE")
    mat <- create_dense_matrix_with_int_dims(100, 200)
    ndarray$write(mat)
    ndarray$close()

    ndarray <- SOMADenseNDArrayOpen(uri)
    ned <- ndarray$non_empty_domain(max_only = TRUE)
    expect_equal(ned, c(99, 199))

    # Test reads out of bounds
    coords <- list(bit64::as.integer64(c(1, 2)), bit64::as.integer64(c(3, 4)))
    expect_no_error(ndarray$read_arrow_table(coords = coords))

    coords <- list(
      bit64::as.integer64(c(101, 202)),
      bit64::as.integer64(c(3, 4))
    )
    expect_error(ndarray$read(coords = coords)$tables()$concat())

    ndarray$close()

    ndarray <- SOMADenseNDArrayOpen(uri, "WRITE")

    # Test resize down
    new_shape <- c(50, 60)
    expect_error(ndarray$resize(new_shape))

    # Test writes out of old bounds
    ndarray <- SOMADenseNDArrayOpen(uri, "WRITE")
    mat <- create_dense_matrix_with_int_dims(300, 400)
    expect_error(ndarray$write(mat))
    ndarray$close()

    # Test resize up, dry run
    new_shape <- c(501, 601)
    expect_no_error(
      reason_string <- ndarray$resize(new_shape, check_only = TRUE)
    )
    expect_equal(reason_string, "")

    # Test resize up, for real
    new_shape <- c(500, 600)
    expect_no_error(ndarray$resize(new_shape))

    # Test writes within new bounds
    ndarray <- SOMADenseNDArrayOpen(uri, "WRITE")
    mat <- create_dense_matrix_with_int_dims(500, 600)
    expect_no_error(ndarray$write(mat))
    ndarray$close()

    ndarray <- SOMADenseNDArrayOpen(uri)
    coords <- list(
      bit64::as.integer64(c(101, 202)),
      bit64::as.integer64(c(3, 4))
    )
    expect_no_condition(x <- ndarray$read_dense_matrix(coords = coords))
    ndarray$close()

    rm(ndarray)
    gc()
  }
})
