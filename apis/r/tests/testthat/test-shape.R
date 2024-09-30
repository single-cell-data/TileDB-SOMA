test_that("SOMADataFrame shape", {
  asch <- create_arrow_schema()

  index_column_name_choices = list(
    "soma_joinid",
    c("soma_joinid", "int_column"),
    c("soma_joinid", "string_column"),
    c("string_column", "int_column"), # duplicate intentional to match `domain_at_create_choices`
    c("string_column", "int_column")
  )

  domain_at_create_choices  = list(
    list(soma_joinid = c(0, 1000)),
    list(soma_joinid = c(0, 1000), int_column = c(-10000, 10000)),
    list(soma_joinid = c(0, 1000), string_column = NULL),
    list(string_column = NULL, int_column = c(-10000, 10000)),
    list(string_column = c("apple", "zebra"), int_column = c(-10000, 10000))
  )

  # Check the test configs themselves to make sure someone (ahem, me)
  # didn't edit one without forgetting to edit the other
  expect_equal(length(index_column_name_choices), length(domain_at_create_choices))

  for (i in seq_along(index_column_name_choices)) {
    index_column_names <- index_column_name_choices[[i]]

    for (use_domain_at_create in c(FALSE, TRUE)) {

      uri <- withr::local_tempdir("soma-dataframe-shape")

      # Create
      if (dir.exists(uri)) unlink(uri, recursive=TRUE)

      domain_for_create <- NULL
      if (use_domain_at_create) {
        domain_for_create <- domain_at_create_choices[[i]]
      }

      sdf <- SOMADataFrameCreate(
        uri,
        asch,
        index_column_names = index_column_names,
        domain = domain_for_create)

      expect_true(sdf$exists())
      expect_true(dir.exists(uri))

      # Write
      tbl0 <- arrow::arrow_table(int_column = 1L:4L,
                                 soma_joinid = 1L:4L,
                                 float_column = 1.1:4.1,
                                 string_column = c("apple", "ball", "cat", "dog"),
                                 schema = asch)

      sdf$write(tbl0)
      sdf$close()

      sdf <- SOMADataFrameOpen(uri)

      # Check shape and maxshape et al.
      if (!.new_shape_feature_flag_is_enabled()) {
        expect_false(sdf$tiledbsoma_has_upgraded_domain())
      } else {
        expect_true(sdf$tiledbsoma_has_upgraded_domain())
      }
      expect_error(sdf$shape(), class = "notYetImplementedError")
      expect_error(sdf$maxshape(), class = "notYetImplementedError")

      # Not implemented this way per
      # https://github.com/single-cell-data/TileDB-SOMA/pull/2953#discussion_r1746125089
      # sjid_shape <- sdf$.maybe_soma_joinid_shape()
      # sjid_maxshape <- sdf$.maybe_soma_joinid_maxshape()
      soma_context <- soma_context()
      sjid_shape <- maybe_soma_joinid_shape(sdf$uri, soma_context)
      sjid_maxshape <- maybe_soma_joinid_maxshape(sdf$uri, soma_context)

      if ("soma_joinid" %in% index_column_names) {
        # More testing to come on
        # https://github.com/single-cell-data/TileDB-SOMA/issues/2407
        expect_false(is.na(sjid_shape))
        expect_false(is.na(sjid_maxshape))
      } else {
        expect_true(is.na(sjid_shape))
        expect_true(is.na(sjid_maxshape))
      }

      # Check has_upgraded_domain
      if (!.new_shape_feature_flag_is_enabled()) {
        expect_false(sdf$tiledbsoma_has_upgraded_domain())
      } else {
        expect_true(sdf$tiledbsoma_has_upgraded_domain())
      }

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

        if (!.new_shape_feature_flag_is_enabled()) {
          # Old behavior
          expect_equal(sjid_dom, sjid_mxd)
        }

        if (!use_domain_at_create) {
          expect_equal(sjid_dom[[1]], 0)
          expect_equal(sjid_mxd[[1]], 0)
          # This is a really big number in the ballpark of 2**63; its exact
          # value is unimportant.
          expect_true(sjid_dom[[2]] > bit64::as.integer64(10000000000))
          expect_true(sjid_mxd[[2]] > bit64::as.integer64(10000000000))
        } else {
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
      }

      if ("int_column" %in% index_column_names) {
        int_dom <- dom[["int_column"]]
        int_mxd <- mxd[["int_column"]]
        int_dfc <- domain_for_create[["int_column"]]

        if (!.new_shape_feature_flag_is_enabled()) {
          # Old behavior
          expect_equal(int_dom, int_mxd)
        }

        if (!use_domain_at_create) {
          expect_true(int_dom[[1]] < -2000000000)
          expect_true(int_dom[[2]] > 2000000000)
        } else {
          expect_equal(int_dom, int_dfc)
        }

        if (!.new_shape_feature_flag_is_enabled()) {
          if (!use_domain_at_create) {
            expect_true(int_mxd[[1]] < -2000000000)
            expect_true(int_mxd[[2]] > 2000000000)
          } else {
            expect_equal(int_mxd, int_dfc)
          }
        } else {
            expect_true(int_mxd[[1]] < -2000000000)
            expect_true(int_mxd[[2]] > 2000000000)
        }
      }

      if ("string_column" %in% index_column_names) {
        str_dom <- dom[["string_column"]]
        str_mxd <- mxd[["string_column"]]
        str_dfc <- domain_for_create[["string_column"]]

        if (!.new_shape_feature_flag_is_enabled()) {
          expect_equal(str_dom, c("", ""))
          expect_equal(str_mxd, c("", ""))

        } else {
          if (!use_domain_at_create) {
            expect_equal(str_dom, c("", ""))
          } else {
            if (is.null(str_dfc)) {
                expect_equal(str_dom, c("", ""))
            } else {
              expect_equal(str_dom, str_dfc)
            }
          }
          expect_equal(str_mxd, c("", ""))
        }
      }

      sdf$close()

      rm(sdf, tbl0)

      gc()
    }
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
    domain = list(soma_joinid = list(bit64::as.integer64(0L), bit64::as.integer64(99L)))
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

test_that("SOMASparseNDArray shape", {
  uri <- withr::local_tempdir("soma-sparse-ndarray-shape")
  asch <- create_arrow_schema()

  element_type_choices = list(arrow::float32(), arrow::int16())
  arg_shape <- c(100, 200)
  for (element_type in element_type_choices) {

    if (dir.exists(uri)) unlink(uri, recursive=TRUE)
    ndarray <- SOMASparseNDArrayCreate(uri, element_type, shape = arg_shape)
    ndarray$close()

    ndarray <- SOMASparseNDArrayOpen(uri)
    expect_equal(ndarray$ndim(), length(arg_shape))

    expect_equal(ndarray$shape(), as.integer64(arg_shape))

    # More generally after current-domain support:
    readback_shape <- ndarray$shape()
    readback_maxshape <- ndarray$maxshape()
    expect_equal(length(readback_shape), length(readback_maxshape))
    if (.new_shape_feature_flag_is_enabled()) {
      expect_true(all(readback_shape < readback_maxshape))
    } else {
      expect_true(all(readback_shape == readback_maxshape))
    }

    ndarray$close()

    # Test write in bounds
    ndarray <- SOMASparseNDArrayOpen(uri, "WRITE")
    soma_dim_0 <- c(2,3)
    soma_dim_1 <- c(4,5)
    soma_data <- c(60, 70)
    sm <- sparseMatrix(i = soma_dim_0, j = soma_dim_1, x = soma_data)
    ndarray$write(sm)
    ndarray$close()

    ndarray <- SOMASparseNDArrayOpen(uri)
    ned <- ndarray$non_empty_domain()
    expect_equal(ned, c(2,4))

    # Test reads out of bounds
    coords <- list(bit64::as.integer64(c(1,2)), bit64::as.integer64(c(3,4)))
    expect_no_error(x <- ndarray$read(coords=coords)$tables()$concat())

    coords <- list(bit64::as.integer64(c(101,202)), bit64::as.integer64(c(3,4)))
    expect_error(x <- ndarray$read(coords=coords)$tables()$concat())

    ndarray$close()

    if (.new_shape_feature_flag_is_enabled()) {
      ndarray <- SOMASparseNDArrayOpen(uri, "WRITE")

      # Test resize down
      new_shape <- c(50, 60)
      expect_error(ndarray$resize(new_shape))

      # Test writes out of old bounds
      soma_dim_0 <- c(200,300)
      soma_dim_1 <- c(400,500)
      soma_data <- c(6000, 7000)
      sm <- sparseMatrix(i = soma_dim_0, j = soma_dim_1, x = soma_data)
      expect_error(ndarray$write(sm))

      # Test resize up
      new_shape <- c(500, 600)
      expect_no_error(ndarray$resize(new_shape))

      # Test writes within new bounds
      soma_dim_0 <- c(200,300)
      soma_dim_1 <- c(400,500)
      soma_data <- c(6000, 7000)
      sm <- sparseMatrix(i = soma_dim_0, j = soma_dim_1, x = soma_data)
      expect_no_error(ndarray$write(sm))
      ndarray$close()

      ndarray <- SOMASparseNDArrayOpen(uri)
      coords <- list(bit64::as.integer64(c(101,202)), bit64::as.integer64(c(3,4)))
      expect_no_error(x <- ndarray$read(coords=coords)$tables()$concat())
      ndarray$close()
    }

    rm(ndarray)
    gc()
  }
})

test_that("SOMADenseNDArray shape", {
  uri <- withr::local_tempdir("soma-dense-ndarray-shape")
  asch <- create_arrow_schema()

  element_type_choices = list(arrow::float32(), arrow::int16())
  arg_shape <- c(100, 200)
  for (element_type in element_type_choices) {

    if (dir.exists(uri)) unlink(uri, recursive=TRUE)
    ndarray <- SOMADenseNDArrayCreate(uri, element_type, shape = arg_shape)
    ndarray$close()

    ndarray <- SOMADenseNDArrayOpen(uri)

    expect_equal(ndarray$ndim(), length(arg_shape))

    expect_equal(ndarray$shape(), as.integer64(arg_shape))

    # More generally after current-domain support:
    readback_shape <- ndarray$shape()
    readback_maxshape <- ndarray$maxshape()
    expect_equal(length(readback_shape), length(readback_maxshape))
    # TODO: Awaiting core support for new shape in dense arrays.
    # https://github.com/single-cell-data/TileDB-SOMA/issues/2955
    #if (.new_shape_feature_flag_is_enabled()) {
    #  expect_true(all(readback_shape < readback_maxshape))
    #} else {
    #  expect_true(all(readback_shape == readback_maxshape))
    #}
    expect_true(all(readback_shape == readback_maxshape))

    ndarray$close()

    # Test write in bounds
    ndarray <- SOMADenseNDArrayOpen(uri, "WRITE")
    mat <- create_dense_matrix_with_int_dims(100, 200)
    ndarray$write(mat)
    ndarray$close()

    ndarray <- SOMADenseNDArrayOpen(uri)
    ned <- ndarray$non_empty_domain()
    expect_equal(ned, c(99, 199))

    # Test reads out of bounds
    coords <- list(bit64::as.integer64(c(1,2)), bit64::as.integer64(c(3,4)))
    expect_no_error(ndarray$read_arrow_table(coords=coords))

    coords <- list(bit64::as.integer64(c(101,202)), bit64::as.integer64(c(3,4)))
    expect_error(ndarray$read(coords=coords)$tables()$concat())

    ndarray$close()

    if (.new_shape_feature_flag_is_enabled()) {
      ndarray <- SOMADenseNDArrayOpen(uri, "WRITE")

      # Test resize down
      new_shape <- c(50, 60)
      expect_error(ndarray$resize(new_shape))

      # Test writes out of old bounds
      ndarray <- SOMADenseNDArrayOpen(uri, "WRITE")
      mat <- create_dense_matrix_with_int_dims(300, 400)
      expect_error(ndarray$write(mat))
      ndarray$close()

      # Test resize up
      new_shape <- c(500, 600)
      # TODO: Awaiting core support for new shape in dense arrays.
      # https://github.com/single-cell-data/TileDB-SOMA/issues/2955
      # expect_no_error(ndarray$resize(new_shape))
      expect_error(ndarray$resize(new_shape))

      # Test writes within new bounds
      ndarray <- SOMADenseNDArrayOpen(uri, "WRITE")
      mat <- create_dense_matrix_with_int_dims(300, 400)
      # TODO: Awaiting core support for new shape in dense arrays.
      # https://github.com/single-cell-data/TileDB-SOMA/issues/2955
      # expect_no_error(ndarray$write(sm))
      expect_error(ndarray$write(sm))
      ndarray$close()

      ndarray <- SOMADenseNDArrayOpen(uri)
      coords <- list(bit64::as.integer64(c(101,202)), bit64::as.integer64(c(3,4)))
      # TODO: Awaiting core support for new shape in dense arrays.
      # https://github.com/single-cell-data/TileDB-SOMA/issues/2955
      # expect_no_error(x <- ndarray$read(coords=coords)$tables()$concat())
      expect_error(x <- ndarray$read(coords=coords)$tables()$concat())
      ndarray$close()
    }

    rm(ndarray)
    gc()
  }
})
