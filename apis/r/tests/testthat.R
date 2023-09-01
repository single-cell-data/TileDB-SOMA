library(testthat)
library(tiledbsoma)

tiledbsoma::show_package_versions()
# more verbose reporting:  test_check("tiledbsoma", reporter=default_reporter())
test_check("tiledbsoma")
