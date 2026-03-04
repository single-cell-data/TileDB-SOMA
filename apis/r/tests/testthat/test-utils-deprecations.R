# .deprecate() ----------------------------------------------------------

test_that(".deprecate() validates inputs", {
  # what must be a length-1 character string
  expect_error(.deprecate(when = "1.0.0", what = 1))
  expect_error(.deprecate(when = "1.0.0", what = c("a", "b")))
  # dots must be empty
  expect_error(.deprecate(when = "1.0.0", what = "foo()", extra = TRUE))
})

test_that(".deprecate() is silent for future deprecations", {
  local_mocked_bindings(.deprecation_stage = function(...) NULL)
  expect_silent(.deprecate(when = "99.0.0", what = "foo()"))
})

test_that(".deprecate() warns during deprecation stage", {
  local_mocked_bindings(.deprecation_stage = function(...) "deprecate")
  expect_warning(
    .deprecate(when = "1.0.0", what = "foo()"),
    class = "lifecycle_warning_deprecated"
  )
})

test_that(".deprecate() errors during defunct stage", {
  local_mocked_bindings(.deprecation_stage = function(...) "defunct")
  expect_error(
    .deprecate(when = "1.0.0", what = "foo()"),
    class = "lifecycle_error_deprecated"
  )
})

# .deprecation_stage() --------------------------------------------------

test_that(".deprecation_stage() validates inputs", {
  # when must be string parseable as a version
  expect_error(.deprecation_stage(1))
  expect_error(.deprecation_stage(""))
  expect_error(
    .deprecation_stage("not-a-version"),
    "'when' must be a valid version"
  )
})

test_that("future release returns NULL", {
  # a version > than all known releases is assumed to be a future release
  # so no deprecation is signaled
  expect_null(.deprecation_stage("99.0.0"))
})

test_that("unknown release errors", {
  # a version that's <= the latest release but not in releases.dcf is invalid
  expect_error(.deprecation_stage("1.99.0"), "Unknown tiledbsoma release")
})

test_that("current version < when returns NULL", {
  # if the current package version hasn't reached the deprecation version yet,
  # no deprecation is signaled
  local_mocked_bindings(.tiledbsoma_deprecation_version = function() {
    package_version("1.16.0")
  })
  expect_null(.deprecation_stage("2.0.0"))
})

test_that("deprecate stage", {
  # same major version and minor diff of 1 triggers deprecation warning
  # (current=2.1.0, when=2.0.0)
  local_mocked_bindings(.tiledbsoma_deprecation_version = function() {
    package_version("2.1.0")
  })
  expect_equal(.deprecation_stage("2.0.0"), "deprecate")
})

test_that("defunct via major version bump", {
  # a major version bump always triggers defunct (current=2.0.0, when=1.17.0)
  local_mocked_bindings(.tiledbsoma_deprecation_version = function() {
    package_version("2.0.0")
  })
  expect_equal(.deprecation_stage("1.17.0"), "defunct")
})

test_that("defunct via minor version + time", {
  # within the same major version, defunct requires both:
  #   - minor version diff >= 2
  #   - >= 12 weeks between release dates
  # 1.15.0 released 2024-12-17, 1.17.0 released 2025-05-21
  # minor diff = 2, ~22 weeks apart (> 12 week threshold)
  local_mocked_bindings(.tiledbsoma_deprecation_version = function() {
    package_version("1.17.0")
  })
  expect_equal(.deprecation_stage("1.15.0"), "defunct")
})
