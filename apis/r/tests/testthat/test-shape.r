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

    expect_false(sdf$has_upgraded_domain())
    expect_error(sdf$shape(), class = "notYetImplementedError")
    expect_error(sdf$maxshape(), class = "notYetImplementedError")

    # Not implemented this way per
    # https://github.com/single-cell-data/TileDB-SOMA/pull/2953#discussion_r1746125089
    # sjid_shape <- sdf$.maybe_soma_joinid_shape()
    # sjid_maxshape <- sdf$.maybe_soma_joinid_maxshape()
    sjid_shape <- maybe_soma_joinid_shape(sdf$uri, config = as.character(tiledb::config(sdf$tiledbsoma_ctx$context())))
    sjid_maxshape <- maybe_soma_joinid_maxshape(sdf$uri, config = as.character(tiledb::config(sdf$tiledbsoma_ctx$context())))

    if (has_soma_joinid_dim) {
      # More testing to come on
      # https://github.com/single-cell-data/TileDB-SOMA/issues/2407
      expect_false(sjid_shape == -1)
      expect_false(sjid_maxshape == -1)
    } else {
      expect_true(sjid_shape == -1)
      expect_true(sjid_maxshape == -1)
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
  shape_choices = list(
    c(100),
    c(100,200),
    c(100,200,300)
  )
  for (element_type in element_type_choices) {
    for (shape in shape_choices) {

      if (dir.exists(uri)) unlink(uri, recursive=TRUE)
      ndarray <- SOMASparseNDArrayCreate(uri, element_type, shape = shape)

      expect_equal(ndarray$ndim(), length(shape))

      expect_equal(ndarray$shape(), as.integer64(shape))

      # More generally after current-domain support:
      readback_shape <- ndarray$shape()
      readback_maxshape <- ndarray$maxshape()
      expect_equal(length(readback_shape), length(readback_maxshape))
      for (i in 1:length(shape)) {
        s <- as.integer(readback_shape[[i]])
        ms <- as.integer(readback_maxshape[[i]])
        expect_true(s <= ms)
      }

      # Resize tests upcoming on
      # https://github.com/single-cell-data/TileDB-SOMA/issues/2407
      ndarray$close()

      rm(ndarray)
      gc()
    }
  }
})
