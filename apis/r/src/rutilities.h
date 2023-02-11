
// Collect helper functions from tiledb-r which it might export via a header

#pragma once

#include <tiledb/tiledb>
#include <spdl.h>

namespace tdbs = tiledbsoma;

// create a single 'comparable' number out of version, minor and patch
#define TileDB_Version(v,m,p)	(((v) * 65536) + ((m) * 256) + (p))

// current build is encoded in TILEDB_VERSION
#define TILEDB_VERSION TileDB_Version(TILEDB_VERSION_MAJOR,TILEDB_VERSION_MINOR,TILEDB_VERSION_PATCH)

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
    std::memcpy(&(num[0]), &(vec[0]), n*sizeof(double));
    return num;
}

// Applies (named list of) vectors of points to the named dimensions
void apply_dim_points(tdbs::SOMAReader* sr,
                      std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>>& name2dim,
                      Rcpp::List lst);

// Applies (named list of) matrices of points to the named dimensions
void apply_dim_ranges(tdbs::SOMAReader* sr,
                      std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>>& name2dim,
                      Rcpp::List lst);
