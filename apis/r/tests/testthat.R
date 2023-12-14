library(testthat)
library(tiledbsoma)

tiledbsoma::show_package_versions()
test_check("tiledbsoma", reporter=default_reporter())
