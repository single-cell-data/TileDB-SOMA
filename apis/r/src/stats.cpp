
#include <Rcpp.h>
#include <tiledb/tiledb>

// [[Rcpp::export(.tiledb_stats_enable)]]
void libtiledb_stats_enable() {
    tiledb::Stats::enable();
}

// [[Rcpp::export(.tiledb_stats_disable)]]
void libtiledb_stats_disable() {
    tiledb::Stats::disable();
}

// [[Rcpp::export(.tiledb_stats_reset)]]
void libtiledb_stats_reset() {
    tiledb::Stats::reset();
}

// [[Rcpp::export(.tiledb_stats_raw_dump)]]
std::string libtiledb_stats_dump() {
    std::string txt;
    tiledb::Stats::raw_dump(&txt);
    return txt;
}
