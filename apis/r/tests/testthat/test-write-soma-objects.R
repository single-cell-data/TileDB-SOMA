
test_that("write_soma.data.frame mechanics", {
  skip_if(!extended_tests())
  skip_if_not_installed('datasets')

  uri <- withr::local_tempdir("write-soma-data-frame")
  collection <- SOMACollectionCreate(uri)

  co2 <- get_data('CO2', package = 'datasets')
  expect_no_condition(sdf <- write_soma(co2, uri = 'co2', soma = collection))
  expect_s3_class(sdf, 'SOMADataFrame')
  expect_true(sdf$exists())
  expect_identical(sdf$uri, file.path(collection$uri, 'co2'))
  expect_identical(sdf$dimnames(), 'soma_joinid')
  expect_identical(sdf$attrnames(), c(names(co2), 'obs_id'))
  expect_true(rlang::is_na(sdf$shape()))
  schema <- sdf$schema()
  expect_s3_class(schema, 'Schema')
  expect_equal(schema$num_fields - 2L, ncol(co2))
  expect_identical(
    setdiff(schema$names, c('soma_joinid', 'obs_id')),
    names(co2)
  )

  collection$close()
})

test_that("write_soma.data.frame enumerations", {
  skip_if(!extended_tests())
  skip_if_not_installed('datasets')

  uri <- withr::local_tempdir("write-soma-data-frame-enumerations")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  co2 <- get_data('CO2', package = 'datasets')
  expect_no_condition(sdf <- write_soma(co2, uri = 'co2', soma = collection))
  expect_s3_class(sdf, 'SOMADataFrame')
  expect_true(sdf$exists())
  expect_identical(sdf$uri, file.path(collection$uri, 'co2'))
  sdf$close()
  expect_no_condition(sdf <- SOMADataFrameOpen(sdf$uri))

  expect_s3_class(schema <- sdf$schema(), 'Schema')
  expect_equal(schema$num_fields - 2L, ncol(co2))
  expect_identical(
    setdiff(schema$names, c('soma_joinid', 'obs_id')),
    names(co2)
  )

  expect_s3_class(tbl <- sdf$read()$concat(), 'Table')
  expect_equal(ncol(tbl), schema$num_fields)
  expect_identical(
    sort(names(tbl)),
    sort(schema$names)
  )

  expect_s3_class(df <- as.data.frame(tbl), 'data.frame')
  for (i in names(co2)) {
    if (is.factor(co2[[i]])) {
      expect_true(is.factor(df[[i]]))
    } else {
      expect_equal(class(df[[i]]), class(co2[[i]]))
    }
  }
})

test_that("write_soma.data.frame no enumerations", {
  skip_if(!extended_tests())
  skip_if_not_installed('datasets')

  uri <- withr::local_tempdir("write-soma-data-frame-factorless")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  co2 <- get_data('CO2', package = 'datasets')
  for (i in names(co2)) {
    if (is.factor(co2[[i]])) {
      co2[[i]] <- as.vector(co2[[i]])
    }
  }

  expect_no_condition(sdf <- write_soma(co2, uri = 'co2', soma = collection))
  expect_s3_class(sdf, 'SOMADataFrame')
  expect_true(sdf$exists())
  expect_identical(sdf$uri, file.path(collection$uri, 'co2'))
  sdf$close()
  expect_no_condition(sdf <- SOMADataFrameOpen(sdf$uri))

  expect_s3_class(schema <- sdf$schema(), 'Schema')
  expect_equal(schema$num_fields - 2L, ncol(co2))
  expect_identical(
    setdiff(schema$names, c('soma_joinid', 'obs_id')),
    names(co2)
  )

  expect_s3_class(tbl <- sdf$read()$concat(), 'Table')
  expect_equal(ncol(tbl), schema$num_fields)
  expect_identical(
    sort(names(tbl)),
    sort(schema$names)
  )

  expect_s3_class(df <- as.data.frame(tbl), 'data.frame')
  for (i in names(co2)) {
    if (is.factor(co2[[i]])) {
      expect_false(is.factor(df[[i]]))
      expect_type(df[[i]], 'character')
    } else {
      expect_equal(class(df[[i]]), class(co2[[i]]))
    }
  }
})

test_that("write_soma.data.frame registration", {
  skip_if(!extended_tests())
  skip_if_not_installed("datasets")

  uri <- withr::local_tempdir("write-soma-data-frame-registration")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  co2 <- get_data("CO2", package = "datasets")

  expect_no_condition(sdf <- write_soma(
    co2,
    uri = "co2",
    soma_parent = collection,
    key = "co2"
  ))
  expect_s3_class(sdf, "SOMADataFrame")
  expect_true(sdf$exists())
  expect_identical(sdf$uri, file.path(collection$uri, "co2"))

  sdf$close()
  collection$close()

  expect_no_condition(collection <- SOMACollectionOpen(collection$uri))
  expect_s3_class(collection, "SOMACollection")
  expect_identical(collection$length(), 1L)
  expect_identical(collection$names(), "co2")
  expect_s3_class(cdf <- collection$get("co2"), "SOMADataFrame")
  expect_s3_class(df <- as.data.frame(cdf$read()$concat()), 'data.frame')
  for (col in names(co2)) {
    expect_identical(df[[col]], co2[[col]])
  }

  # Registration assertions
  expect_error(write_soma(co2, "uri", soma_parent = collection, key = TRUE))
  expect_error(write_soma(co2, "uri", soma_parent = collection, key = 1L))
  expect_error(write_soma(co2, "uri", soma_parent = collection, key = c("a", "b")))
  expect_error(write_soma(co2, "uri", soma_parent = NULL, key = "co2"))

})

test_that("write_soma dense matrix mechanics", {
  skip_if(!extended_tests())
  skip_if_not_installed('datasets')

  uri <- withr::local_tempdir("write-soma-dense-matrix")
  collection <- SOMACollectionCreate(uri)

  state77 <- get(x = 'state.x77', envir = getNamespace('datasets'))
  expect_no_condition(dmat <- write_soma(
    state77,
    uri = 'state77',
    soma = collection,
    sparse = FALSE
  ))
  expect_s3_class(dmat, 'SOMADenseNDArray')
  expect_true(dmat$exists())
  expect_identical(dmat$uri, file.path(collection$uri, 'state77'))
  expect_equal(dmat$ndim(), 2L)
  expect_identical(dmat$dimnames(), paste0('soma_dim_', c(0L, 1L)))
  expect_identical(dmat$attrnames(), 'soma_data')
  expect_equal(dmat$shape(), dim(state77))
  # Test transposition
  expect_no_condition(tmat <- write_soma(
    state77,
    uri = 'state77t',
    soma = collection,
    sparse = FALSE,
    transpose = TRUE
  ))
  expect_s3_class(tmat, 'SOMADenseNDArray')
  expect_true(tmat$exists())
  expect_identical(tmat$uri, file.path(collection$uri, 'state77t'))
  expect_equal(tmat$ndim(), 2L)
  expect_identical(tmat$dimnames(), paste0('soma_dim_', c(0L, 1L)))
  expect_identical(tmat$attrnames(), 'soma_data')
  expect_equal(tmat$shape(), rev(dim(state77)))
  # Error if given sparse matrix and ask for dense
  knex <- get_data('KNex', package = 'Matrix')$mm
  expect_error(write_soma(knex, uri = 'knex', soma = collection, sparse = FALSE))
  # Work on dgeMatrices
  expect_no_condition(emat <- write_soma(
    as(knex, 'unpackedMatrix'),
    uri = 'knexd',
    soma = collection,
    sparse = FALSE
  ))
  expect_s3_class(emat, 'SOMADenseNDArray')
  expect_true(emat$exists())
  expect_identical(emat$uri, file.path(collection$uri, 'knexd'))
  expect_equal(emat$ndim(), 2L)
  expect_identical(emat$dimnames(), paste0('soma_dim_', c(0L, 1L)))
  expect_identical(emat$attrnames(), 'soma_data')
  expect_equal(emat$shape(), dim(knex))

  collection$close()
})

test_that("write_soma dense matrix registration", {
  skip_if(!extended_tests())
  skip_if_not_installed("datasets")

  uri <- withr::local_tempdir("write-soma-dense-matrix-registration")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  state77 <- get(x = 'state.x77', envir = getNamespace('datasets'))

  expect_no_condition(dmat <- write_soma(
    state77,
    uri = "state77",
    soma_parent = collection,
    sparse = FALSE,
    key = "state77"
  ))
  expect_s3_class(dmat, "SOMADenseNDArray")
  expect_true(dmat$exists())
  expect_identical(dmat$uri, file.path(collection$uri, "state77"))

  dmat$close()
  collection$close()

  expect_no_condition(collection <- SOMACollectionOpen(collection$uri))
  expect_s3_class(collection, "SOMACollection")
  expect_identical(collection$length(), 1L)
  expect_identical(collection$names(), "state77")
  expect_s3_class(cmat <- collection$get("state77"), "SOMADenseNDArray")
  expect_type(mat <- cmat$read_dense_matrix(), 'double')
  expect_true(is.matrix(mat))
  expect_identical(mat, unname(state77))

  # Registration assertions
  expect_error(write_soma(state77, "uri", soma_parent = collection, key = TRUE))
  expect_error(write_soma(state77, "uri", soma_parent = collection, key = 1L))
  expect_error(write_soma(state77, "uri", soma_parent = collection, key = c("a", "b")))
  expect_error(write_soma(state77, "uri", soma_parent = NULL, key = "state77"))

})

test_that("write_soma sparse matrix mechanics", {
  skip_if(!extended_tests())
  uri <- withr::local_tempdir("write-soma-sparse-matrix")
  collection <- SOMACollectionCreate(uri)
  knex <- get_data('KNex', package = 'Matrix')$mm
  expect_no_condition(smat <- write_soma(knex, uri = 'knex', soma = collection))
  expect_s3_class(smat, 'SOMASparseNDArray')
  expect_true(smat$exists())
  expect_identical(smat$uri, file.path(collection$uri, 'knex'))
  expect_equal(smat$ndim(), 2L)
  expect_identical(smat$dimnames(), paste0('soma_dim_', c(0L, 1L)))
  expect_identical(smat$attrnames(), 'soma_data')
  expect_equal(smat$shape(), dim(knex))
  # Test transposition
  expect_no_condition(tmat <- write_soma(
    knex,
    uri = 'knext',
    soma = collection,
    transpose = TRUE
  ))
  expect_s3_class(tmat, 'SOMASparseNDArray')
  expect_true(tmat$exists())
  expect_identical(tmat$uri, file.path(collection$uri, 'knext'))
  expect_equal(tmat$ndim(), 2L)
  expect_identical(tmat$dimnames(), paste0('soma_dim_', c(0L, 1L)))
  expect_identical(tmat$attrnames(), 'soma_data')
  expect_equal(tmat$shape(), rev(dim(knex)))
  # Try a dense matrix
  skip_if_not_installed('datasets')
  state77 <- get(x = 'state.x77', envir = getNamespace('datasets'))
  expect_no_condition(cmat <- write_soma(state77, 'state77s', soma = collection))
  expect_s3_class(cmat, 'SOMASparseNDArray')
  expect_true(cmat$exists())
  expect_identical(cmat$uri, file.path(collection$uri, 'state77s'))
  expect_equal(cmat$ndim(), 2L)
  expect_identical(cmat$dimnames(), paste0('soma_dim_', c(0L, 1L)))
  expect_identical(cmat$attrnames(), 'soma_data')
  expect_equal(cmat$shape(), dim(state77))

  collection$close()
})

test_that("write_soma sparse matrix registration", {
  skip_if(!extended_tests())
  skip_if_not_installed("datasets")

  uri <- withr::local_tempdir("write-sparse-dense-matrix-registration")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  knex <- get_data('KNex', package = 'Matrix')$mm

  expect_no_condition(smat <- write_soma(
    knex,
    uri = "knex",
    soma_parent = collection,
    key = "knex"
  ))
  expect_s3_class(smat, "SOMASparseNDArray")
  expect_true(smat$exists())
  expect_identical(smat$uri, file.path(collection$uri, "knex"))

  smat$close()
  collection$close()

  expect_no_condition(collection <- SOMACollectionOpen(collection$uri))
  expect_s3_class(collection, "SOMACollection")
  expect_identical(collection$length(), 1L)
  expect_identical(collection$names(), "knex")
  expect_s3_class(cmat <- collection$get("knex"), "SOMASparseNDArray")
  expect_s4_class(mat <- cmat$read()$sparse_matrix()$concat(), "dgTMatrix")
  expect_identical(as.matrix(mat), as.matrix(unname(knex)))

  # Registration assertions
  expect_error(write_soma(knex, "uri", soma_parent = collection, key = TRUE))
  expect_error(write_soma(knex, "uri", soma_parent = collection, key = 1L))
  expect_error(write_soma(knex, "uri", soma_parent = collection, key = c("a", "b")))
  expect_error(write_soma(knex, "uri", soma_parent = NULL, key = "knex"))

})

test_that("write_soma.character mechanics", {
  skip_if(!extended_tests())
  skip_if_not_installed("datasets")

  uri <- withr::local_tempdir("write-soma-character")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  mtcars <- get_data("mtcars", package = "datasets")
  cars <- row.names(mtcars)

  expect_no_condition(sdf <- write_soma(cars, "cars", soma_parent = collection))
  expect_s3_class(sdf, "SOMADataFrame")
  expect_true(sdf$exists())
  expect_length(sdf$attrnames(), 1L)
  expect_identical(sdf$attrnames(), "values")

  hint <- "soma_uns_outgest_hint"
  expect_true(hint %in% names(sdf$get_metadata()))
  expect_identical(sdf$get_metadata(hint), uns_hint("1d")[[hint]])

  sdf$close()
  expect_no_condition(sdf <- SOMADataFrameOpen(sdf$uri))
  expect_s3_class(tbl <- sdf$read()$concat(), "Table")
  expect_identical(nrow(tbl), length(cars))
  expect_identical(tbl$values$as_vector(), cars)
})

test_that("write_soma.character scalar", {
  skip_if(!extended_tests())
  skip_if_not_installed("datasets")

  uri <- withr::local_tempdir("write-soma-character-scalar")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  mtcars <- get_data("mtcars", package = "datasets")
  cars <- paste(row.names(mtcars), collapse = ",")

  expect_no_condition(sdf <- write_soma(cars, "cars", soma_parent = collection))
  expect_s3_class(sdf, "SOMADataFrame")
  expect_true(sdf$exists())
  expect_length(sdf$attrnames(), 1L)
  expect_identical(sdf$attrnames(), "values")

  hint <- "soma_uns_outgest_hint"
  expect_true(hint %in% names(sdf$get_metadata()))
  expect_identical(sdf$get_metadata(hint), uns_hint("1d")[[hint]])

  sdf$close()
  expect_no_condition(sdf <- SOMADataFrameOpen(sdf$uri))
  expect_s3_class(tbl <- sdf$read()$concat(), "Table")
  expect_identical(nrow(tbl), 1L)
  expect_identical(tbl$soma_joinid$as_vector(), 0L)
  expect_identical(sdf.cars <- tbl$values$as_vector(), cars)
  expect_length(unlist(strsplit(sdf.cars, ",")), nrow(mtcars))
})
