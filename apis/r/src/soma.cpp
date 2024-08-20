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

    auto soup = tdbs::SOMAObject::open(uri, OpenMode::read, sctx);
    auto tpstr = soup->type();
    if (!tpstr.has_value()) {
        Rcpp::stop("No object type value for URI '%s'", uri);
    }
    return tpstr.value();
}

// [[Rcpp::export]]
std::string get_tiledb_object_type(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // shared pointer to TileDB Context from SOMAContext
    std::shared_ptr<tiledb::Context> ctx = sctx->tiledb_ctx();

    auto objtype = tiledb::Object::object(*ctx, uri).type();
    switch (objtype) {
        case tdbs::Object::Type::Array:
            return std::string("ARRAY");
            break;
        case tdbs::Object::Type::Group:
            return std::string("GROUP");
            break;
        case tdbs::Object::Type::Invalid:
            return std::string("INVALID");
            break;
        default:
            Rcpp::stop("Inadmissable object type ('%d') for URI '%s'", (int)objtype, uri);
            break;
    }
    Rcpp::stop("No object type value for URI '%s'", uri);
}
