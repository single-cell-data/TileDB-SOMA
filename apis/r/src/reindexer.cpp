#include <Rcpp/Lightest>        // for R interface to C++
#include <RcppInt64>            // for fromInteger64

#include <tiledbsoma/tiledbsoma>
#if TILEDB_VERSION_MAJOR == 2 && TILEDB_VERSION_MINOR >= 4
#include <tiledb/tiledb_experimental>
#endif
#include <tiledbsoma/reindexer/reindexer.h>

namespace tdbs = tiledbsoma;

#include "xptr-utils.h"         // xptr taggging utilitie


// [[Rcpp::export]]
Rcpp::XPtr<tdbs::IntIndexer> reindex_create() {
    auto p = new tdbs::IntIndexer();
    return make_xptr<tdbs::IntIndexer>(p);
}

// [[Rcpp::export]]
Rcpp::XPtr<tdbs::IntIndexer> reindex_map(Rcpp::XPtr<tdbs::IntIndexer> idx,
                                         const Rcpp::NumericVector nvec) {
    check_xptr_tag<tdbs::IntIndexer>(idx);
    const std::vector<int64_t> vec = Rcpp::fromInteger64(nvec);
    idx->map_locations(vec);
    return idx;
}

// [[Rcpp::export]]
Rcpp::NumericVector reindex_lookup(Rcpp::XPtr<tdbs::IntIndexer> idx,
                                   const Rcpp::NumericVector kvec) {
    check_xptr_tag<tdbs::IntIndexer>(idx);
    const std::vector<int64_t> keys = Rcpp::fromInteger64(kvec);
    int sz = keys.size();
    std::vector<int64_t> res(sz);
    idx->lookup(keys, res);
    return Rcpp::toInteger64(res);
}
