// we currently get deprecation warnings by default which are noisy
#ifndef TILEDB_NO_API_DEPRECATION_WARNINGS
#define TILEDB_NO_API_DEPRECATION_WARNINGS
#endif

#include <Rcpp.h>
#include <tiledbsoma/tiledbsoma>

//' TileDB Statistics interface
//'
//' The functions `tiledbsoma_stats_enable`, `tiledbsoma_stats_disable`, `tiledbsoma_stats_reset`
//' and `tiledbsoma_stats_dump` expose the TileDB Core functionality for performance measurements
//' and statistics.  The first three just turn on, off or reset, the fourth returns a JSON string.
//' For convenience the function `tiledbsoma_stats_show` displays the information on the console.
//'
//' @export
// [[Rcpp::export]]
void tiledbsoma_stats_enable() {
    tiledbsoma::stats::enable();
}

//' @rdname tiledbsoma_stats_enable
//' @export
// [[Rcpp::export]]
void tiledbsoma_stats_disable() {
    tiledbsoma::stats::disable();
}

//' @rdname tiledbsoma_stats_enable
//' @export
// [[Rcpp::export]]
void tiledbsoma_stats_reset() {
    tiledbsoma::stats::reset();
}

//' @rdname tiledbsoma_stats_enable
//' @export
// [[Rcpp::export]]
std::string tiledbsoma_stats_dump() {
    return tiledbsoma::stats::dump();
}
