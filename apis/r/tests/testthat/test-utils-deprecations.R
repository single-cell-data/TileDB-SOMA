test_that("input type validation", {
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
