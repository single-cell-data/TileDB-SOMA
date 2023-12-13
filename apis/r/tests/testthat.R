if (Sys.info()[["sysname"]] == "Linux") {
    Sys.setenv("TILEDB_SM_COMPUTE_CONCURRENCY_LEVEL"=1,
               "TILEDB_SM_IO_CONCURRENCY_LEVEL"=1,
               "OMP_THREAD_LIMIT"=1)
}

library(testthat)
library(tiledbsoma)

tiledbsoma::show_package_versions()
test_check("tiledbsoma", reporter=default_reporter())
