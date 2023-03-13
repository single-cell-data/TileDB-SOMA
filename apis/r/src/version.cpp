
#include <Rcpp.h>
#include <tiledbsoma/tiledbsoma>

//' libtiledbsoma version information
//'
//' Returns a string with version information for libtiledbsoma and the linked TileDB Embedded library.
// [[Rcpp::export]]
std::string libtiledbsoma_version() {
    return tiledbsoma::version::as_string();
}

//' TileDB Embedded Version interface
//'
//' This gets the version of the TileDB Embedded library that is currently in use.
// [[Rcpp::export]]
Rcpp::IntegerVector tiledb_embedded_version() {
    std::tuple<int, int, int> triple = tiledbsoma::version::embedded_version_triple();
    return Rcpp::IntegerVector::create(std::get<0>(triple), std::get<1>(triple), std::get<2>(triple));
}
