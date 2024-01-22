##spdl::set_level('warn')

test_that("Write Assay mechanics", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))

  uri <- withr::local_tempdir("write-assay")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close())

  rna <- get_data('pbmc_small', package = 'SeuratObject')[['RNA']]
  expect_no_condition(ms <- write_soma(rna, soma_parent = collection))
  on.exit(ms$close())
  expect_s3_class(ms, 'SOMAMeasurement')
  expect_true(ms$exists())

  expect_identical(ms$uri, file.path(collection$uri, 'rna'))
  expect_identical(ms$names(), c('X', 'var'))
  expect_s3_class(ms$var, 'SOMADataFrame')
  expect_identical(setdiff(ms$var$attrnames(), 'var_id'), names(rna[[]]))
  expect_s3_class(ms$X, 'SOMACollection')
  layers <- c(counts = 'counts', data = 'data', scale.data = 'scale_data')
  expect_identical(ms$X$names(), unname(layers))
  for (i in seq_along(layers)) {
    expect_equal(
      ms$X$get(layers[i])$shape(),
      rev(dim(rna)),
      info = layers[i]
    )
  }

  # Test no feature-level meta data
  rna2 <- rna
  for (i in names(rna2[[]])) {
    rna2[[i]] <- NULL
  }
  expect_no_condition(ms2 <- write_soma(
    rna2,
    uri = 'rna-no-md',
    soma_parent = collection
  ))
  expect_s3_class(ms2, 'SOMAMeasurement')
  expect_true(ms2$exists())
  expect_identical(ms2$uri, file.path(collection$uri, 'rna-no-md'))
  expect_identical(ms2$names(), c('X', 'var'))
  expect_s3_class(ms2$var, 'SOMADataFrame')
  expect_identical(ms2$var$attrnames(), 'var_id')
  ms2$close()

  # Test no counts
  rna3 <- SeuratObject::SetAssayData(rna, 'counts', new('matrix'))
  expect_no_condition(ms3 <- write_soma(
    rna3,
    uri = 'rna-no-counts',
    soma_parent = collection
  ))
  expect_s3_class(ms3, 'SOMAMeasurement')
  expect_true(ms3$exists())
  expect_identical(ms3$uri, file.path(collection$uri, 'rna-no-counts'))
  expect_identical(ms3$names(), c('X', 'var'))
  expect_s3_class(ms3$X, 'SOMACollection')
  lyrs <- layers[c('data', 'scale.data')]
  expect_identical(ms3$X$names(), unname(lyrs))
  for (i in seq_along(lyrs)) {
    expect_equal(
      ms3$X$get(lyrs[i])$shape(),
      rev(dim(rna3)),
      info = lyrs[i]
    )
  }
  ms3$close()

  # Test no scale.data
  rna4 <- SeuratObject::SetAssayData(rna, 'scale.data', new('matrix'))
  expect_no_condition(ms4 <- write_soma(
    rna4,
    uri = 'rna-no-scale',
    soma_parent = collection
  ))
  expect_s3_class(ms4, 'SOMAMeasurement')
  expect_true(ms4$exists())
  expect_identical(ms4$uri, file.path(collection$uri, 'rna-no-scale'))
  expect_identical(ms4$names(), c('X', 'var'))
  expect_s3_class(ms4$X, 'SOMACollection')
  lyrs <- layers[c('counts', 'data')]
  expect_identical(ms4$X$names(), unname(lyrs))
  for (i in seq_along(lyrs)) {
    expect_equal(
      ms4$X$get(lyrs[i])$shape(),
      rev(dim(rna4)),
      info = lyrs[i]
    )
  }
  ms4$close()

  # Test no counts or scale.data
  rna5 <- SeuratObject::SetAssayData(rna3, 'scale.data', new('matrix'))
  expect_no_condition(ms5 <- write_soma(
    rna5,
    uri = 'rna-no-counts-scale',
    soma_parent = collection)
  )
  on.exit(ms5$close())
  expect_s3_class(ms5, 'SOMAMeasurement')
  expect_true(ms5$exists())
  expect_identical(ms5$uri, file.path(collection$uri, 'rna-no-counts-scale'))
  expect_identical(ms5$names(), c('X', 'var'))
  expect_s3_class(ms5$X, 'SOMACollection')
  lyrs <- layers[c('counts', 'data')]
  expect_identical(ms5$X$names(), 'data')
  expect_equal(ms5$X$get('data')$shape(), rev(dim(rna5)))

  # Test assertions
  expect_error(write_soma(rna, uri = TRUE, soma_parent = collection))
  expect_error(write_soma(rna, uri = c('dir', 'rna'), soma_parent = collection))
  expect_error(write_soma(
    rna,
    soma_parent = SOMADataFrameCreate(uri = file.path(uri, 'data-frame'))
  ))

  # Verify data slot isn't ingested when it's identical to counts
  rna6 <- SeuratObject::CreateAssayObject(
    counts = SeuratObject::GetAssayData(rna, "counts")
  )
  expect_identical(
    SeuratObject::GetAssayData(rna6, "counts"),
    SeuratObject::GetAssayData(rna6, "data")
  )
  expect_no_condition(ms6 <- write_soma(
    rna6,
    uri = "rna-identical-counts-data",
    soma_parent = collection
  ))
  expect_equal(ms6$X$names(), "counts")
})

test_that("Write DimReduc mechanics", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))

  uri <- withr::local_tempdir("write-reduction")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close())
  pbmc_small <- get_data('pbmc_small', package = 'SeuratObject')
  pbmc_small_rna <- pbmc_small[['RNA']]
  pbmc_small_pca <- pbmc_small[['pca']]
  pbmc_small_tsne <- pbmc_small[['tsne']]

  # Test writing PCA
  ms_pca <- write_soma(pbmc_small_rna, uri = 'rna-pca', soma_parent = collection)
  on.exit(ms_pca$close())
  fidx <- match(rownames(SeuratObject::Loadings(pbmc_small_pca)), rownames(pbmc_small_rna))
  expect_no_condition(write_soma(
    pbmc_small_pca,
    soma_parent = ms_pca,
    fidx = fidx,
    nfeatures = nrow(pbmc_small_rna)
  ))
  expect_true(all(ms_pca$names() %in% c('X', 'var', 'obsm', 'varm')))
  expect_identical(ms_pca$obsm$names(), 'X_pca')
  expect_s3_class(spca <- ms_pca$obsm$get('X_pca'), 'SOMASparseNDArray')
  expect_equal(spca$shape(), dim(pbmc_small_pca))
  expect_identical(ms_pca$varm$names(), 'PCs')
  expect_s3_class(sldgs <- ms_pca$varm$get('PCs'), 'SOMASparseNDArray')
  expect_equal(sldgs$shape(), c(nrow(pbmc_small_rna), ncol(pbmc_small_pca)))

  # Test writing tSNE
  ms_tsne <- write_soma(pbmc_small_rna, uri = 'rna-tsne', soma_parent = collection)
  on.exit(ms_tsne$close())
  expect_no_condition(write_soma(pbmc_small_tsne, soma_parent = ms_tsne))
  expect_true(all(ms_tsne$names() %in% c('X', 'var', 'obsm', 'varm')))
  expect_identical(ms_tsne$obsm$names(), 'X_tsne')
  expect_s3_class(stsne <- ms_tsne$obsm$get('X_tsne'), 'SOMASparseNDArray')
  expect_equal(stsne$shape(), dim(pbmc_small_tsne))

  # Test writing both PCA and tSNE
  ms <- write_soma(pbmc_small_rna, soma_parent = collection)
  expect_no_condition(ms_pca2 <- write_soma(
    pbmc_small_pca,
    soma_parent = ms,
    fidx = fidx,
    nfeatures = nrow(pbmc_small_rna)
  ))
  on.exit(ms_pca2$close())
  expect_no_condition(ms_rna <- write_soma(pbmc_small_rna, soma_parent = ms))
  on.exit(ms_rna$close())
  expect_true(all(ms$names() %in% c('X', 'var', 'obsm', 'varm')))
  expect_true(all(ms$obsm$names() %in% paste0('X_', c('pca', 'tsne'))))
  expect_identical(ms$varm$names(), 'PCs')

  # Test assertions
  expect_error(write_soma(pbmc_small_pca, uri = 'X_pca', soma_parent = ms_tsne))
  expect_error(write_soma(pbmc_small_pca, soma_parent = collection))
  expect_true(ms_tsne$is_open())
  expect_warning(ms_pca3 <- write_soma(pbmc_small_pca, soma_parent = ms_tsne))
  expect_error(ms_tsne$varm)
})

test_that("Write Graph mechanics", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))

  uri <- withr::local_tempdir("write-graph")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close())

  pbmc_small <- get_data('pbmc_small', package = 'SeuratObject')
  pbmc_small_rna <- pbmc_small[['RNA']]
  graph <- pbmc_small[['RNA_snn']]

  ms <- write_soma(pbmc_small_rna, soma_parent = collection)
  on.exit(ms$close())
  expect_no_condition(write_soma(graph, uri = 'rna-snn', soma_parent = ms))
  expect_true(all(ms$names() %in% c('X', 'var', 'obsp')))
  expect_identical(ms$obsp$names(), 'rna-snn')
  expect_s3_class(sgrph <- ms$obsp$get('rna-snn'), 'SOMASparseNDArray')
  expect_equal(sgrph$shape(), dim(graph))

  # Test assertions
  expect_error(write_soma(grph, collection = soma_parent))

  expect_identical(ms$uri, file.path(collection$uri, 'rna'))
  expect_identical(ms$names(), c('X', 'var', 'obsp'))
  expect_s3_class(ms$var, 'SOMADataFrame')
  expect_identical(setdiff(ms$var$attrnames(), 'var_id'), names(pbmc_small_rna[[]]))
  expect_s3_class(ms$X, 'SOMACollection')
  layers <- c(counts = 'counts', data = 'data', scale.data = 'scale_data')
  expect_identical(ms$X$names(), unname(layers))
  for (i in seq_along(layers)) {
    expect_equal(
      ms$X$get(layers[i])$shape(),
      rev(dim(pbmc_small_rna)),
      info = layers[i]
    )
  }

  # Test no feature-level meta data
  rna2 <- pbmc_small_rna
  for (i in names(rna2[[]])) {
    rna2[[i]] <- NULL
  }
  expect_no_condition(ms2 <- write_soma(
    rna2,
    uri = 'rna-no-md',
    soma_parent = collection
  ))
  on.exit(ms2$close())
  expect_s3_class(ms2, 'SOMAMeasurement')
  expect_true(ms2$exists())
  expect_identical(ms2$uri, file.path(collection$uri, 'rna-no-md'))
  expect_identical(ms2$names(), c('X', 'var'))
  expect_s3_class(ms2$var, 'SOMADataFrame')
  expect_identical(ms2$var$attrnames(), 'var_id')

  # Test no counts
  rna3 <- SeuratObject::SetAssayData(pbmc_small_rna, 'counts', new('matrix'))
  expect_no_condition(ms3 <- write_soma(
    rna3,
    uri = 'rna-no-counts',
    soma_parent = collection
  ))
  on.exit(ms3$close())
  expect_s3_class(ms3, 'SOMAMeasurement')
  expect_true(ms3$exists())
  expect_identical(ms3$uri, file.path(collection$uri, 'rna-no-counts'))
  expect_identical(ms3$names(), c('X', 'var'))
  expect_s3_class(ms3$X, 'SOMACollection')
  lyrs <- layers[c('data', 'scale.data')]
  expect_identical(ms3$X$names(), unname(lyrs))
  for (i in seq_along(lyrs)) {
    expect_equal(
      ms3$X$get(lyrs[i])$shape(),
      rev(dim(rna3)),
      info = lyrs[i]
    )
  }

  # Test no scale.data
  rna4 <- SeuratObject::SetAssayData(pbmc_small_rna, 'scale.data', new('matrix'))
  expect_no_condition(ms4 <- write_soma(
    rna4,
    uri = 'rna-no-scale',
    soma_parent = collection
  ))
  on.exit(ms4$close())
  expect_s3_class(ms4, 'SOMAMeasurement')
  expect_true(ms4$exists())
  expect_identical(ms4$uri, file.path(collection$uri, 'rna-no-scale'))
  expect_identical(ms4$names(), c('X', 'var'))
  expect_s3_class(ms4$X, 'SOMACollection')
  lyrs <- layers[c('counts', 'data')]
  expect_identical(ms4$X$names(), unname(lyrs))
  for (i in seq_along(lyrs)) {
    expect_equal(
      ms4$X$get(lyrs[i])$shape(),
      rev(dim(rna4)),
      info = lyrs[i]
    )
  }

  # Test no counts or scale.data
  rna5 <- SeuratObject::SetAssayData(rna3, 'scale.data', new('matrix'))
  expect_no_condition(ms5 <- write_soma(
    rna5,
    uri = 'rna-no-counts-scale',
    soma_parent = collection)
  )
  on.exit(ms5$close())
  expect_s3_class(ms5, 'SOMAMeasurement')
  expect_true(ms5$exists())
  expect_identical(ms5$uri, file.path(collection$uri, 'rna-no-counts-scale'))
  expect_identical(ms5$names(), c('X', 'var'))
  expect_s3_class(ms5$X, 'SOMACollection')
  lyrs <- layers[c('counts', 'data')]
  expect_identical(ms5$X$names(), 'data')
  expect_equal(ms5$X$get('data')$shape(), rev(dim(rna5)))

  # Test assertions
  expect_error(write_soma(pbmc_small_rna, uri = TRUE, soma_parent = collection))
  expect_error(write_soma(pbmc_small_rna, uri = c('dir', 'rna'), soma_parent = collection))
  expect_error(write_soma(
    pbmc_small_rna,
    soma_parent = SOMADataFrameCreate(uri = file.path(uri, 'data-frame'))
  ))
})

test_that("Write DimReduc mechanics", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("write-reduction-2")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close())
  pbmc_small <- get_data('pbmc_small', package = 'SeuratObject')
  pbmc_small_rna <- pbmc_small[['RNA']]
  pbmc_small_pca <- pbmc_small[['pca']]
  pbmc_small_tsne <- pbmc_small[['tsne']]

  # Test writing PCA
  mspca <- write_soma(pbmc_small_rna, uri = 'rna-pca', soma_parent = collection)
  on.exit(mspca$close())
  fidx <- match(rownames(SeuratObject::Loadings(pbmc_small_pca)), rownames(pbmc_small_rna))
  expect_no_condition(write_soma(
    pbmc_small_pca,
    soma_parent = mspca,
    fidx = fidx,
    nfeatures = nrow(pbmc_small_rna)
  ))
  expect_true(all(mspca$names() %in% c('X', 'var', 'obsm', 'varm')))
  expect_identical(mspca$obsm$names(), 'X_pca')
  expect_s3_class(spca <- mspca$obsm$get('X_pca'), 'SOMASparseNDArray')
  expect_equal(spca$shape(), dim(pbmc_small_pca))
  expect_identical(mspca$varm$names(), 'PCs')
  expect_s3_class(sldgs <- mspca$varm$get('PCs'), 'SOMASparseNDArray')
  expect_equal(sldgs$shape(), c(nrow(pbmc_small_rna), ncol(pbmc_small_pca)))

  # Test writing tSNE
  mstsne <- write_soma(pbmc_small_rna, uri = 'rna-tsne', soma_parent = collection)
  on.exit(mstsne$close())
  expect_no_condition(write_soma(pbmc_small_tsne, soma_parent = mstsne))
  expect_true(all(mstsne$names() %in% c('X', 'var', 'obsm', 'varm')))
  expect_identical(mstsne$obsm$names(), 'X_tsne')
  expect_s3_class(stsne <- mstsne$obsm$get('X_tsne'), 'SOMASparseNDArray')
  expect_equal(stsne$shape(), dim(pbmc_small_tsne))

  # Test writing both PCA and tSNE
  ms <- write_soma(pbmc_small_rna, soma_parent = collection)
  on.exit(ms$close())
  expect_no_condition(write_soma(
    pbmc_small_pca,
    soma_parent = ms,
    fidx = fidx,
    nfeatures = nrow(pbmc_small_rna)
  ))
  expect_no_condition(write_soma(pbmc_small_rna, soma_parent = ms))
  expect_true(all(ms$names() %in% c('X', 'var', 'obsm', 'varm')))
  expect_true(all(ms$obsm$names() %in% paste0('X_', c('pca', 'tsne'))))
  expect_identical(ms$varm$names(), 'PCs')

  # Test assertions
  expect_error(write_soma(pbmc_small_pca, uri = 'X_pca', soma_parent = mstsne))
  expect_error(write_soma(pbmc_small_pca, soma_parent = collection))
  expect_warning(write_soma(pbmc_small_pca, soma_parent = mstsne))
  expect_error(mstsne$varm)
})

test_that("Write Graph mechanics", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("write-graph")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close())
  pbmc_small <- get_data('pbmc_small', package = 'SeuratObject')
  pbmc_small_rna <- pbmc_small[['RNA']]
  graph <- pbmc_small[['RNA_snn']]
  ms <- write_soma(pbmc_small_rna, soma_parent = collection)
  on.exit(ms$close())
  expect_no_condition(write_soma(graph, uri = 'rna-snn', soma_parent = ms))
  expect_true(all(ms$names() %in% c('X', 'var', 'obsp')))
  expect_identical(ms$obsp$names(), 'rna-snn')
  expect_s3_class(sgrph <- ms$obsp$get('rna-snn'), 'SOMASparseNDArray')
  expect_equal(sgrph$shape(), dim(graph))

  # Test assertions
  expect_error(write_soma(grph, collection = soma_parent))
})

test_that("Write Seurat mechanics", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  pbmc_small <- get_data('pbmc_small', package = 'SeuratObject')
  uri <- withr::local_tempdir(SeuratObject::Project(pbmc_small))

  expect_no_condition(uri <- write_soma(pbmc_small, uri))

  expect_type(uri, 'character')
  expect_true(grepl(
    paste0('^', SeuratObject::Project(pbmc_small)),
    basename(uri)
  ))

  expect_no_condition(experiment <- SOMAExperimentOpen(uri))
  on.exit(experiment$close())

  expect_no_error(experiment$ms)

  expect_equal(experiment$mode(), "READ")
  expect_s3_class(experiment, 'SOMAExperiment')
  expect_true(grepl(
    paste0('^', SeuratObject::Project(pbmc_small)),
    basename(experiment$uri)
  ))

  expect_identical(experiment$ms$names(), 'RNA')
  expect_s3_class(ms <- experiment$ms$get('RNA'), 'SOMAMeasurement')

  expect_identical(
    lapply(list(ms$X$names()), sort),
    lapply(list(c('counts', 'data', 'scale_data')), sort)
  )
  expect_identical(
    lapply(list(ms$obsm$names()), sort),
    lapply(list(c('X_pca', 'X_tsne')), sort)
  )
  expect_identical(ms$varm$names(), 'PCs')
  expect_identical(ms$obsp$names(), 'RNA_snn')
  expect_error(ms$varp)
  expect_identical(
    setdiff(experiment$obs$attrnames(), 'obs_id'),
    names(pbmc_small[[]])
  )
})
