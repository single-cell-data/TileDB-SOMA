
#include <Rcpp.h>
#include <tiledbsoma/tiledbsoma>

//' TileDB Embedded Version interface
//'
//' This gets the version of the TileDB Embedded library that is currently in use.
std::string tiledb_embedded_version() {
    return tiledbsoma::version::as_string();
}
