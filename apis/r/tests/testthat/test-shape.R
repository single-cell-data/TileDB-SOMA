test_that("SOMADataFrame shape", {
  #uri <- withr::local_tempdir("soma-dataframe-shape")
  uri <- "/tmp/fooze"
  asch <- create_arrow_schema()

  index_column_name_choices = list(
    "soma_joinid",
    c("soma_joinid", "foo"),
    c("soma_joinid", "baz"),
    c("baz", "foo")
  )

  for (index_column_names in index_column_name_choices) {
    has_soma_joinid_dim <- "soma_joinid" %in% index_column_names

    if (dir.exists(uri)) unlink(uri, recursive=TRUE)

    sdf <- SOMADataFrameCreate(uri, asch, index_column_names = index_column_names)
    expect_true(sdf$exists())
    expect_true(dir.exists(uri))

    tbl0 <- arrow::arrow_table(foo = 1L:4L,
                               soma_joinid = 1L:4L,
                               bar = 1.1:4.1,
                               baz = c("apple", "ball", "cat", "dog"),
                               schema = asch)

    sdf$write(tbl0)
    sdf$close()

    sdf <- SOMADataFrameOpen(uri)

    if (.new_shape_feature_flag_is_enabled()) {
      expect_true(sdf$tiledbsoma_has_upgraded_domain())
    } else {
      expect_false(sdf$tiledbsoma_has_upgraded_domain())
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

    if (has_soma_joinid_dim) {
      # More testing to come on
      # https://github.com/single-cell-data/TileDB-SOMA/issues/2407
      expect_false(is.na(sjid_shape))
      expect_false(is.na(sjid_maxshape))
    } else {
      expect_true(is.na(sjid_shape))
      expect_true(is.na(sjid_maxshape))
    }

    sdf$close()

    rm(sdf, tbl0)

    gc()
  }
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
