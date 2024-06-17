test_that("Load SeuratComand mechanics", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  skip_if_not_installed('jsonlite')

  pbmc_small <- get_data('pbmc_small', package = 'SeuratObject')
  uri <- write_soma(pbmc_small, uri = tempfile(pattern='load-seurat-command'))

  expect_no_condition(exp <- SOMAExperimentOpen(uri))
  on.exit(exp$close(), add = TRUE)

  expect_s3_class(uns <- exp$get('uns'), 'SOMACollection')
  expect_in('seurat_commands', uns$names())
  expect_s3_class(logs <- uns$get('seurat_commands'), 'SOMACollection')
  expect_identical(sort(logs$names()), sort(SeuratObject::Command(pbmc_small)))

  expect_length(
    cmds <- .load_seurat_command(uns, SeuratObject::DefaultAssay(pbmc_small)),
    length(SeuratObject::Command(pbmc_small))
  )
  expect_identical(sort(names(cmds)), sort(SeuratObject::Command(pbmc_small)))

  for (cmd in names(cmds)) {
    expect_s4_class(log <- cmds[[cmd]], 'SeuratCommand')
    objlog <- pbmc_small[[cmd]]
    if (!is.null(SeuratObject::DefaultAssay(log))) {
      expect_identical(
        SeuratObject::DefaultAssay(log),
        SeuratObject::DefaultAssay(pbmc_small)
      )
    }
    for (i in setdiff(slotNames(log), 'params')) {
      info <- paste("Differing slot", sQuote(i), "for command", sQuote(cmd))
      switch(
        EXPR = i,
        time.stamp = expect_equal(
          methods::slot(log, i),
          methods::slot(objlog, i),
          info = info
        ),
        call.string = {
          objstring <- paste(trimws(methods::slot(objlog, i)), collapse = ' ')
          expect_identical(methods::slot(log, i), objstring, info = info)
        },
        expect_identical(methods::slot(log, i), methods::slot(objlog, i), info = info)
      )
    }
    params <- methods::slot(log, 'params')
    for (param in names(params)) {
      info <- paste("Differing parameter", sQuote(param), "for command", sQuote(cmd))
      objparam <- do.call(`$`, list(objlog, param))
      switch(
        typeof(params[[param]]),
        character = {
          objstring <- paste(trimws(objparam), collapse = ' ')
          expect_identical(params[[param]], objstring, info = info)
        },
        double = expect_equal(params[[param]], objparam, info = info),
        expect_identical(params[[param]], objparam, info = info)
      )
    }
  }

  no.assay <- Filter(
    function(x) is.null(SeuratObject::DefaultAssay(x)),
    methods::slot(pbmc_small, 'commands')
  )
  expect_length(no.cmds <- .load_seurat_command(uns, 'no-assay'), length(no.assay))
  expect_identical(sort(names(no.cmds)), sort(names(no.assay)))

  # Test assertions
  expect_error(.load_seurat_command(list(), SeuratObject::DefaultAssay(pbmc_small)))
  expect_error(.load_seurat_command(exp$ms, SeuratObject::DefaultAssay(pbmc_small)))
  expect_error(.load_seurat_command(uns, 1:3), regexp = "^'ms_names' must be a character vector with no empty strings$")
  expect_error(.load_seurat_command(uns, TRUE), regexp = "^'ms_names' must be a character vector with no empty strings$")
  expect_error(.load_seurat_command(uns, NULL), regexp = "^'ms_names' must be a character vector with no empty strings$")
  expect_error(.load_seurat_command(uns, ''), regexp = "^'ms_names' must be a character vector with no empty strings$")
})

test_that("Loading SeuratCommands works from experiment queries", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  skip_if_not_installed('jsonlite')


  pbmc_small <- get_data('pbmc_small', package = 'SeuratObject')
  uri <- write_soma(pbmc_small, uri = tempfile(pattern='seurat-command-query'))

  expect_no_condition(exp <- SOMAExperimentOpen(uri))
  on.exit(exp$close(), add = TRUE)

  expect_s3_class(
    query <- SOMAExperimentAxisQuery$new(exp, SeuratObject::DefaultAssay(pbmc_small)),
    'SOMAExperimentAxisQuery'
  )

  expect_s4_class(
    obj <- suppressWarnings(query$to_seurat(X_layers = c('data' = 'data'))),
    'Seurat'
  )
  expect_identical(
    sort(SeuratObject::Command(obj)),
    sort(SeuratObject::Command(pbmc_small))
  )
})

test_that("Load SeuratCommand with missing commands", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  skip_if_not_installed('jsonlite')

  pbmc_small <- get_data('pbmc_small', package = 'SeuratObject')
  slot(pbmc_small, "commands") <- list()
  expect_true(validObject(pbmc_small))
  uri <- write_soma(pbmc_small, uri = tempfile(pattern='missing-commands'))

  expect_no_condition(exp <- SOMAExperimentOpen(uri))
  on.exit(exp$close(), add = TRUE)
  expect_true('uns' %in% exp$names())
  expect_s3_class(uns <- exp$get('uns'), 'SOMACollection')
  expect_false('seurat_commands' %in% uns$names())

  expect_s3_class(
    query <- SOMAExperimentAxisQuery$new(exp, SeuratObject::DefaultAssay(pbmc_small)),
    'SOMAExperimentAxisQuery'
  )

  expect_no_condition(obj <- query$to_seurat(X_layers = c('data' = 'data')))

  withr::with_options(
    list(verbose = TRUE),
    expect_warning(
      query$to_seurat(X_layers = c('data' = 'data')),
      regexp = "^Cannot find a SOMACollection with command logs in 'uns'$"
    )
  )

  expect_s4_class(obj, 'Seurat')
  expect_true(validObject(obj))
  expect_length(SeuratObject::Command(obj), 0L)
})
