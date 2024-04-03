
// Collect helper functions from tiledb-r which it might export via a header

#pragma once

#include <spdl.h>
#include <tiledb/tiledb>

namespace tdbs = tiledbsoma;

// create a single 'comparable' number out of version, minor and patch
#define TileDB_Version(v, m, p) (((v)*65536) + ((m)*256) + (p))

// current build is encoded in TILEDB_VERSION
#define TILEDB_VERSION TileDB_Version(TILEDB_VERSION_MAJOR, \
                                      TILEDB_VERSION_MINOR, \
                                      TILEDB_VERSION_PATCH)

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

inline ResultOrder get_tdb_result_order(std::string result_order){
	std::map<std::string, ResultOrder> result_order_map{
		{"auto", ResultOrder::automatic},
		{"row-major", ResultOrder::rowmajor},
		{"column-major", ResultOrder::colmajor}
	};
	return result_order_map[result_order];
}

struct ContextWrapper {
    //ContextWrapper(std::shared_ptr<tiledb::Context> ctx_ptr_) : ctxptr(std::move(ctx_ptr_)) {}
    ContextWrapper(std::shared_ptr<tiledb::Context> ctx_ptr_) : ctxptr(ctx_ptr_) {}
    std::shared_ptr<tiledb::Context> ctxptr;
};
typedef struct ContextWrapper ctx_wrap_t;

inline void exitIfError(const ArrowErrorCode ec, const std::string& msg) {
    if (ec != NANOARROW_OK) Rcpp::stop(msg);
}

// Attaches a schema to an array external pointer. The nanoarrow R package
// attempts to do this whenever possible to avoid misinterpreting arrays.
inline void array_xptr_set_schema(SEXP array_xptr, SEXP schema_xptr) {
    R_SetExternalPtrTag(array_xptr, schema_xptr);
}
