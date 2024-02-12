test_that("Write Assay mechanics", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))

  uri <- withr::local_tempdir("write-assay")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)

  rna <- get_data('pbmc_small', package = 'SeuratObject')[['RNA']]
  expect_no_condition(ms <- write_soma(rna, soma_parent = collection))
  on.exit(ms$close(), add = TRUE, after = FALSE)
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
      expected = rev(dim(rna)),
      info = layers[i]
    )
  }
  # Test no feature-level meta data
  rna2 <- rna
  for (i in names(rna2[[]])) {
    rna2[[i]] <- NULL
  }
  expect_no_condition(ms2 <- write_soma(rna2, uri = 'rna-no-md', soma_parent = collection))
  on.exit(ms2$close(), add = TRUE, after = FALSE)
  expect_s3_class(ms2, 'SOMAMeasurement')
  expect_true(ms2$exists())
  expect_identical(ms2$uri, file.path(collection$uri, 'rna-no-md'))
  expect_identical(ms2$names(), c('X', 'var'))
  expect_s3_class(ms2$var, 'SOMADataFrame')
  expect_identical(ms2$var$attrnames(), 'var_id')

  # Test no counts
  rna3 <- SeuratObject::SetAssayData(rna, 'counts', new('matrix'))
  expect_no_condition(ms3 <- write_soma(rna3, uri = 'rna-no-counts', soma_parent = collection))
  on.exit(ms3$close(), add = TRUE, after = FALSE)
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
      expected = rev(dim(rna3)),
      info = lyrs[i]
    )
  }

  # Test no scale.data
  rna4 <- SeuratObject::SetAssayData(rna, 'scale.data', new('matrix'))
  expect_no_condition(ms4 <- write_soma(rna4, uri = 'rna-no-scale', soma_parent = collection))
  on.exit(ms4$close(), add = TRUE, after = FALSE)
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
      expected = rev(dim(rna4)),
      info = lyrs[i]
    )
  }

  # Test no counts or scale.data
  rna5 <- SeuratObject::SetAssayData(rna3, 'scale.data', new('matrix'))
  expect_no_condition(ms5 <- write_soma(rna5, uri = 'rna-no-counts-scale', soma_parent = collection))
  on.exit(ms5$close(), add = TRUE, after = FALSE)
  expect_s3_class(ms5, 'SOMAMeasurement')
  expect_true(ms5$exists())
  expect_identical(ms5$uri, file.path(collection$uri, 'rna-no-counts-scale'))
  expect_identical(ms5$names(), c('X', 'var'))
  expect_s3_class(ms5$X, 'SOMACollection')
  lyrs <- layers[c('counts', 'data')]
  expect_identical(ms5$X$names(), 'data')
  expect_equal(ms5$X$get('data')$shape(), rev(dim(rna5)))

  # Verify data slot isn't ingested when it's identical to counts
  rna6 <- SeuratObject::CreateAssayObject(
    counts = SeuratObject::GetAssayData(rna, "counts")
  )
  expect_identical(
    SeuratObject::GetAssayData(rna6, "counts"),
    SeuratObject::GetAssayData(rna6, "data")
  )
  expect_no_condition(ms6 <- write_soma(rna6, uri = "rna-identical-counts-data", soma_parent = collection))
  on.exit(ms6$close(), add = TRUE, after = FALSE)
  expect_equal(ms6$X$names(), "counts")

  # Test assertions
  expect_error(write_soma(rna, uri = TRUE, soma_parent = collection))
  expect_error(write_soma(rna, uri = c('dir', 'rna'), soma_parent = collection))
  expect_error(write_soma(
    rna,
    soma_parent = SOMADataFrameCreate(uri = file.path(uri, 'data-frame'))
  ))
})

test_that("Write DimReduc mechanics", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))

  uri <- withr::local_tempdir("write-reduction")
  collection <- SOMACollectionCreate(uri)
  on.exit(collection$close(), add = TRUE)
  pbmc_small <- get_data('pbmc_small', package = 'SeuratObject')
  pbmc_small_rna <- pbmc_small[['RNA']]
  pbmc_small_pca <- pbmc_small[['pca']]
  pbmc_small_tsne <- pbmc_small[['tsne']]

  # Test writing PCA
  ms_pca <- write_soma(pbmc_small_rna, uri = 'rna-pca', soma_parent = collection)
  on.exit(ms_pca$close(), add = TRUE, after = FALSE)
  fidx <- match(rownames(SeuratObject::Loadings(pbmc_small_pca)), rownames(pbmc_small_rna))
  expect_no_condition(write_soma(
    pbmc_small_pca,
    soma_parent = ms_pca,
    fidx = fidx,
    nfeatures = nrow(pbmc_small_rna)
  ))
  expect_identical(sort(ms_pca$names()), sort(c('X', 'var', 'obsm', 'varm')))
  expect_identical(ms_pca$obsm$names(), 'X_pca')
  expect_s3_class(spca <- ms_pca$obsm$get('X_pca'), 'SOMASparseNDArray')
  expect_equal(spca$shape(), dim(pbmc_small_pca))
  expect_identical(ms_pca$varm$names(), 'PCs')
  expect_s3_class(sldgs <- ms_pca$varm$get('PCs'), 'SOMASparseNDArray')
  expect_equal(sldgs$shape(), c(nrow(pbmc_small_rna), ncol(pbmc_small_pca)))

  # Test writing tSNE
  ms_tsne <- write_soma(pbmc_small_rna, uri = 'rna-tsne', soma_parent = collection)
  on.exit(ms_tsne$close(), add = TRUE, after = FALSE)
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
  on.exit(ms_pca2$close(), add = TRUE, after = FALSE)
  expect_no_condition(write_soma(pbmc_small_tsne, soma_parent = ms))
  expect_identical(sort(ms$names()), sort(c('X', 'var', 'obsm', 'varm')))
  expect_identical(sort(ms$obsm$names()), sort(paste0('X_', c('pca', 'tsne'))))
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
  on.exit(collection$close(), add = TRUE)

  pbmc_small <- get_data('pbmc_small', package = 'SeuratObject')
  pbmc_small_rna <- pbmc_small[['RNA']]
  graph <- pbmc_small[['RNA_snn']]

  ms <- write_soma(pbmc_small_rna, soma_parent = collection)
  on.exit(ms$close(), add = TRUE, after = FALSE)
  expect_no_condition(write_soma(graph, uri = 'rna-snn', soma_parent = ms))
  expect_identical(sort(ms$names()), sort(c('X', 'var', 'obsp')))
  expect_identical(ms$obsp$names(), 'rna-snn')
  expect_s3_class(sgrph <- ms$obsp$get('rna-snn'), 'SOMASparseNDArray')
  expect_equal(sgrph$shape(), dim(graph))

  # Test assertions
  expect_error(write_soma(graph, collection = soma_parent))
})

test_that("Write SeuratCommand mechanics", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  skip_if_not_installed('jsonlite')

  uri <- withr::local_tempdir('write-command-log')
  uns <- SOMACollectionCreate(uri)
  on.exit(uns$close, add = TRUE)

  pbmc_small <- get_data('pbmc_small', package = 'SeuratObject')
  for (cmd in SeuratObject::Command(pbmc_small)) {
    cmdlog <- pbmc_small[[cmd]]
    cmdlist <- as.list(cmdlog)
    # Test dumping the command log to SOMA
    expect_no_condition(write_soma(cmdlog, uri = cmd, soma_parent = uns), )
    expect_s3_class(cmdgrp <- uns$get('seurat_commands'), 'SOMACollection')

    expect_s3_class(cmddf <- cmdgrp$get(cmd), 'SOMADataFrame')
    expect_identical(cmddf$mode(), 'CLOSED')

    expect_no_condition(cmddf <- cmddf$open('READ', internal_use_only = 'allowed_use'))
    on.exit(cmddf$close(), add = TRUE, after = FALSE)

    # Test qualities of the SOMADataFrame
    expect_identical(cmddf$attrnames(), 'values')
    expect_identical(sort(cmddf$colnames()), sort(c('soma_joinid', 'values')))
    expect_identical(basename(cmddf$uri), cmd)
    expect_equal(cmddf$ndim(), 1L)

    # Test reading the SOMADataFrame
    expect_s3_class(tbl <- cmddf$read()$concat(), 'Table')
    expect_equal(dim(tbl), c(1L, 2L))
    expect_identical(colnames(tbl), cmddf$colnames())
    expect_s3_class(df <- as.data.frame(tbl), 'data.frame')
    expect_type(df$values, 'character')

    # Test decoding the JSON-encoded command log
    expect_type(vals <- jsonlite::fromJSON(df$values), 'list')
    # Test slots of the command log
    for (slot in setdiff(methods::slotNames(cmdlog), 'params')) {
      cmdslot <- methods::slot(cmdlog, slot)
      cmdslot <- if (is.null(cmdslot)) {
        cmdslot
      } else if (inherits(cmdslot, 'POSIXt')) {
        cmdslot <- as.character(jsonlite::toJSON(
          sapply(
            unclass(as.POSIXlt(cmdslot)),
            .encode_as_char,
            simplify = FALSE,
            USE.NAMES = TRUE
          ),
          auto_unbox = TRUE
        ))
      } else if (is.character(cmdslot)) {
        paste(trimws(cmdslot), collapse = ' ')
      } else {
        as.character(cmdslot)
      }
      expect_identical(vals[[slot]], cmdslot)
    }
    # Test encoded parameters
    expect_length(params <- vals[names(cmdlist)], length(cmdlist))
    expect_identical(sort(names(params)), sort(names(cmdlist)))
    for (param in names(params)) {
      if (is.character(cmdlist[[param]])) {
        expect_identical(params[[param]], cmdlist[[param]])
      } else if (is.double(cmdlist[[param]])) {
        # Doubles are encoded as hexadecimal
        expect_identical(params[[param]], sprintf('%a', cmdlist[[param]]))
      } else {
        expect_equivalent(params[[param]], cmdlist[[param]])
      }
    }
  }
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

  expect_s3_class(experiment, 'SOMAExperiment')
  expect_equal(experiment$mode(), "READ")
  expect_true(grepl(
    paste0('^', SeuratObject::Project(pbmc_small)),
    basename(experiment$uri)
  ))

  expect_no_error(experiment$ms)
  expect_identical(experiment$ms$names(), 'RNA')
  expect_s3_class(ms <- experiment$ms$get('RNA'), 'SOMAMeasurement')

  expect_identical(sort(ms$X$names()), sort(c('counts', 'data', 'scale_data')))
  expect_identical(sort(ms$obsm$names()), sort(c('X_pca', 'X_tsne')))
  expect_identical(ms$varm$names(), 'PCs')
  expect_identical(ms$obsp$names(), 'RNA_snn')
  expect_error(ms$varp)
  expect_identical(
    setdiff(experiment$obs$attrnames(), 'obs_id'),
    names(pbmc_small[[]])
  )

  # Test assertions
  expect_error(write_soma(pbmc_small, TRUE))
  expect_error(write_soma(pbmc_small, 1))
  expect_error(write_soma(pbmc_small, ''))
})
