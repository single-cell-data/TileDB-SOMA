
#include <Rcpp.h>
#include <tiledb/tiledb>

//' TileDB Statistics interface
//'
//' The functions `tiledbsoma_stats_enable`, `tiledbsoma_stats_disable`, `tiledbsoma_stats_reset`
//' and `tiledbsoma_stats_dump` expose the TileDB Core functionality for performance measurements
//' and statistics.  The first three just turn on, off or reset, the fourth returns a JSON string.
//' For convenience the function `tiledbsoma_stats_show` displays the information on the console
//'
//' @export
// [[Rcpp::export]]
void stats_enable() {
    tiledb::Stats::enable();
}

//' @rdname stats_enable
//' @export
// [[Rcpp::export]]
void stats_disable() {
    tiledb::Stats::disable();
}

//' @rdname stats_enable
//' @export
// [[Rcpp::export]]
void stats_reset() {
    tiledb::Stats::reset();
}

//' @rdname stats_enable
//' @export
// [[Rcpp::export]]
std::string stats_dump() {
    std::string txt;
    tiledb::Stats::raw_dump(&txt);
    return txt;
}
