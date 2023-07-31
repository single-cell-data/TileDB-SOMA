# Focuses on https://github.com/single-cell-data/TileDB-SOMA/pull/1524

test_that("membership-caching", {

  uri <- withr::local_tempdir('membership-caching')

  # make exp, open for write
  expect_no_condition(exp <- SOMAExperimentCreate(uri))

  # write exp$ms
  expect_no_condition(exp$ms <- SOMACollectionCreate(file.path(uri, "ms"))$close())
  expect_no_condition(exp$close())

  # exp is open for write
  exp <- SOMAExperimentOpen(uri, "WRITE")
  on.exit(exp$close())
  expect_true(exp$mode() == "WRITE")

  # add exp$ms$get("one")
  ms <- exp$ms
  expect_no_condition(ms$set(SOMAMeasurementCreate(file.path(uri, "ms", "one"))$close(), 'one'))

  # add exp$ms$get("one")$obsm
  expect_true(exp$ms$length() == 1)
  expect_no_condition(one <- exp$ms$get("one"))
  expect_no_condition(one$set(SOMACollectionCreate(file.path(uri, "ms", "one", "obsm"))$close(), 'obsm'))
  expect_true(exp$ms$get("one")$length() == 1)

})
