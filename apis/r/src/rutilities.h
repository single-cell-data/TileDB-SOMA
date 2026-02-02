
// Collect helper functions from tiledb-r which it might export via a header

#pragma once

#include <Rcpp.h>  // for R interface to C++
#include <sstream>
#include <tiledb/tiledb>
#include "tiledbsoma_types.h"
namespace tdbs = tiledbsoma;

// create a single 'comparable' number out of version, minor and patch
#define TileDB_Version(v, m, p) (((v) * 65536) + ((m) * 256) + (p))

// current build is encoded in TILEDB_VERSION
#define TILEDB_VERSION TileDB_Version(TILEDB_VERSION_MAJOR, TILEDB_VERSION_MINOR, TILEDB_VERSION_PATCH)

// Applies (named list of) vectors of points to the named dimensions
void apply_dim_points(
    tdbs::ManagedQuery* mq,
    std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>>& name2dim,
    Rcpp::List lst);

// Applies (named list of) matrices of points to the named dimensions
void apply_dim_ranges(
    tdbs::ManagedQuery* mq,
    std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>>& name2dim,
    Rcpp::List lst);

// Convert R config vector to map<string,string> suitable for SOMAArray
inline std::map<std::string, std::string> config_vector_to_map(Rcpp::Nullable<Rcpp::CharacterVector> config) {
    std::map<std::string, std::string> platform_config;

    if (!config.isNull()) {
        Rcpp::CharacterVector confvec(config.get());
        Rcpp::CharacterVector namesvec = confvec.attr("names");  // extract names from named R vector
        size_t n = confvec.length();
        for (size_t i = 0; i < n; i++) {
            platform_config.emplace(std::make_pair(std::string(namesvec[i]), std::string(confvec[i])));
            std::stringstream ss;
            ss << "[config_vector_to_map] adding '" << std::string(namesvec[i]) << "' = '" << std::string(confvec[i])
               << "'";
            tdbs::common::logging::LOG_TRACE(ss.str());
        }
    }

    return platform_config;
}

inline ResultOrder get_tdb_result_order(std::string result_order) {
    std::map<std::string, ResultOrder> result_order_map{
        {"auto", ResultOrder::automatic},
        {"row-major", ResultOrder::rowmajor},
        {"column-major", ResultOrder::colmajor}};
    return result_order_map[result_order];
}

inline void exitIfError(const ArrowErrorCode ec, const std::string& msg) {
    if (ec != NANOARROW_OK)
        Rcpp::stop(msg);
}

// Attaches a schema to an array external pointer. The nanoarrow R package
// attempts to do this whenever possible to avoid misinterpreting arrays.
inline void array_xptr_set_schema(SEXP array_xptr, SEXP schema_xptr) {
    R_SetExternalPtrTag(array_xptr, schema_xptr);
}

// make a TimestampRange from a DatetimeVector (of size two)
std::optional<tdbs::TimestampRange> makeTimestampRange(Rcpp::Nullable<Rcpp::DatetimeVector> tsvec);

std::vector<int64_t> i64_from_rcpp_numeric(const Rcpp::NumericVector& input);

// Code reuse for non_empty_domain, domain, and maxdomain all of which:
// * call a libtiledbsoma function
// * obtain an ArrowTable
// * need to map that to an R list of lo/hi pairs
SEXP convert_domainish(const tdbs::common::arrow::ArrowTable& arrow_table);

// Maps e.g. "int8" and "float32" to "c" and "f".
std::string remap_arrow_type_code_r_to_c(std::string input);

// Maps tiledb layouts to string identifiers
const char* _tiledb_layout_to_string(tiledb_layout_t layout);

// Get options for a TileDB filter
Rcpp::List _get_filter_options(Rcpp::XPtr<tiledb::Filter> filter);

// Get domain from a TileDB dimension
SEXP _get_dim_domain(Rcpp::XPtr<tiledb::Dimension> dim);

// Get tiling from a TileDB dimension
SEXP _get_dim_tile(Rcpp::XPtr<tiledb::Dimension> dim);
