extended_tests <- function() {
  ## check if at CI, if so extended test
  ## could add if pre-release number ie 1.4.3.1 instead of 1.4.3
  ci_set <- Sys.getenv("CI", "") != ""
  ## check for macOS
  macos <- Sys.info()["sysname"] == "Darwin"
  ## check for possible override of 'force' or 'Force'
  ci_override <- tolower(Sys.getenv("CI", "")) == "force"
  ## run extended tests if CI is set, and if either not macOS or 'force' has been set
  ## (ie setting 'force' will enable on macOS too)
  ci_set && (!macos || ci_override)
}

covr_tests <- function() {
  ## check if coverage is flagged
  ## could add if pre-release number ie 1.4.3.1 instead of 1.4.3
  Sys.getenv("COVR", "") != ""
}
