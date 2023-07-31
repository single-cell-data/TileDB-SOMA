
// Collect helper functions from tiledb-r which it might export via a header

#pragma once

#include <spdl.h>
#include <tiledb/tiledb>

namespace tdbs = tiledbsoma;

// create a single 'comparable' number out of version, minor and patch
#define TileDB_Version(v, m, p) (((v)*65536) + ((m)*256) + (p))

// current build is encoded in TILEDB_VERSION
#define TILEDB_VERSION \
    TileDB_Version(    \
        TILEDB_VERSION_MAJOR, TILEDB_VERSION_MINOR, TILEDB_VERSION_PATCH)

// Create a integer64 object
//
// Integer64 is an S3 class. Integers in R are 32-bits. To handle C++
// signed 64-bit integers (int64_t), the full bits may be stored using double '
// as an intermediary which then can be coereced to Integer64.
// For more on this see e.g.
// https://gallery.rcpp.org/articles/creating-integer64-and-nanotime-vectors/
//
inline Rcpp::NumericVector makeInteger64(const std::vector<int64_t>& vec) {
    size_t n = vec.size();

    Rcpp::NumericVector num(n);
    std::memcpy(&(num[0]), vec.data(), n * sizeof(double));

    num.attr("class") = "integer64";
    return (num);
}

// Convert to a scalar int64_t
//
inline int64_t makeScalarInteger64(const double val) {
    int64_t newval;
    memcpy(&newval, &val, sizeof(double));
    return newval;
}

// Create a int64_t vector from a NumericVector
//
inline std::vector<int64_t> getInt64Vector(Rcpp::NumericVector vec) {
    size_t n = vec.size();
    std::vector<int64_t> num(n);
    std::memcpy(&(num[0]), &(vec[0]), n * sizeof(double));
    return num;
}

// Applies (named list of) vectors of points to the named dimensions
void apply_dim_points(
    tdbs::SOMAArray* sr,
    std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>>&
        name2dim,
    Rcpp::List lst);

// Applies (named list of) matrices of points to the named dimensions
void apply_dim_ranges(
    tdbs::SOMAArray* sr,
    std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>>&
        name2dim,
    Rcpp::List lst);

// Convert R config vector to map<string,string> suitable for SOMAArray
inline std::map<std::string, std::string> config_vector_to_map(Rcpp::Nullable<Rcpp::CharacterVector> config) {
    std::map<std::string, std::string> platform_config;

    if (!config.isNull()) {
        Rcpp::CharacterVector confvec(config.get());
        Rcpp::CharacterVector namesvec = confvec.attr("names"); // extract names from named R vector
        size_t n = confvec.length();
        for (size_t i = 0; i<n; i++) {
            platform_config.emplace(std::make_pair(std::string(namesvec[i]), std::string(confvec[i])));
            spdl::trace("[config_vector_to_map] adding '{}' = '{}'", std::string(namesvec[i]), std::string(confvec[i]));
        }
    }

    return platform_config;
}
