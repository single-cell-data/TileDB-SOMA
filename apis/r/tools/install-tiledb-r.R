#!/usr/bin/env Rscript

options(bspm.version.check = TRUE)

# Find tiledb-r constraints
deps <- read.dcf('DESCRIPTION', fields = 'Imports')[1, ]
deps <- trimws(strsplit(gsub(pattern = ',', replacement = '\n', x = deps), split = '\n')[[1L]])
deps <- Filter(nzchar, deps)
deps <- Filter(function(x) grepl('^tiledb', x), deps)
deps <- gsub(pattern = '[[:space:]]', replacement = '', x = deps)
deps <- gsub(pattern = '.*\\(|\\)', replacement = '', x = deps)
const <- data.frame(comp = gsub(pattern = '[[:digit:]\\.]*', replacement = '', x = deps))
const$vers <- gsub(pattern = paste(const$comp, collapse = '|'), replacement = '', x = deps)
(const)

# Find correct version of tiledb-r
db <- utils::available.packages(filters = c("R_version", "OS_type", "subarch"))
idx <- which(rownames(db) == 'tiledb')
valid <- vapply(
  X = idx,
  FUN = function(i) {
    v <- db[i, 'Version']
    res <- vector(mode = 'logical', length = nrow(const))
    for (i in seq_along(res)) {
      res[i] <- do.call(const$comp[i], args = list(v, const$vers[i]))
    }
    return(all(res))
  },
  FUN.VALUE = logical(1L)
)
if (!any(valid)) {
  stop("No valid versions of tiledb-r found")
}
if (!all(valid)) {
  db <- db[-idx[!valid], , drop = FALSE]
}

# Install upstream deps
(ups <- Reduce(
  f = union,
  x = apply(
    X = db[db[, 'Package'] == 'tiledb', , drop = FALSE],
    MARGIN = 1L,
    FUN = \(x) tools::package_dependencies("tiledb", db = t(as.matrix(x)), recursive = TRUE)$tiledb,
    simplify = FALSE
  )
))
utils::install.packages(intersect(ups, rownames(db)))

# Install correct version of tiledb-r
ctb <- contrib.url(getOption("repos"))
names(ctb) <- getOption("repos")
dbr <- db[idx[valid], 'Repository']
(repos <- names(ctb[ctb %in% dbr]))

# BSPM doesn't respect `repos`
# Check to see if any repo w/ valid tiledb-r is CRAN
# If not, turn of BSPM
cran <- getOption("repos")['CRAN']
cran[is.na(cran)] <- ""
if (requireNamespace("bspm", quietly = TRUE) && !cran %in% repos) {
  bspm::disable()
}
utils::install.packages("tiledb", repos = repos)
