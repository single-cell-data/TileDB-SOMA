#include <Rcpp/Lighter>                         // for R interface to C++
#include <nanoarrow/r.h>                        // for C/C++ interface to Arrow (via header exported from the R package)
#include <nanoarrow/nanoarrow.hpp>              // for C/C++ interface to Arrow (vendored)
#include <RcppInt64>                            // for fromInteger64

#include <tiledbsoma/tiledbsoma>
#include <tiledbsoma/reindexer/reindexer.h>

#include "rutilities.h"         				// local declarations
#include "xptr-utils.h"         				// xptr taggging utilities

namespace tdbs = tiledbsoma;

// [[Rcpp::export]]
std::string get_soma_object_type(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // shared pointer to TileDB Context from SOMAContext
    std::shared_ptr<tiledb::Context> ctx = sctx->tiledb_ctx();

    auto tiledb_type = tdbs::Object::object(*ctx, uri).type();
    switch (tiledb_type) {
        case tdbs::Object::Type::Array:
            return std::string("SOMAArray");
            break;
        case tdbs::Object::Type::Group:
            return std::string("SOMAGroup");
            break;
        default:
            Rcpp::stop("Inadmissable object type for URI '%s'", uri);
            break;
    }
}
