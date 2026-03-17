#include <Rcpp/Lighter>  // for R interface to C++
#include <sstream>

#include <nanoarrow/r.h>            // for C/C++ interface to Arrow (via header exported from the R package)
#include <RcppInt64>                // for fromInteger64
#include <nanoarrow/nanoarrow.hpp>  // for C/C++ interface to Arrow (vendored)

#include <tiledbsoma/reindexer/reindexer.h>
#include <tiledbsoma/tiledbsoma>

#include "metadata.h"
#include "rutilities.h"  // local declarations
#include "xptr-utils.h"  // xptr taggging utilities

namespace tdbs = tiledbsoma;

// [[Rcpp::export]]
void c_group_create(
    std::string& uri,
    std::string& type,
    Rcpp::XPtr<somactx_wrap_t> ctxxp,
    Rcpp::Nullable<Rcpp::DatetimeVector> timestamp = R_NilValue) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // shared pointer to TileDB Context from SOMAContext
    // not needed here:  std::shared_ptr<tiledb::Context> ctx =
    // sctx->tiledb_ctx();

    // optional timestamp range
    std::optional<tdbs::TimestampRange> tsrng = makeTimestampRange(timestamp);
    std::stringstream ss;
    ss << "[c_group_create] uri " << uri;
    if (timestamp.isNotNull()) {
        Rcpp::DatetimeVector v(timestamp);
        ss << " ts (" << v[0] << ", " << v[1] << ")";
    }
    tdbs::common::logging::LOG_DEBUG(ss.str());

    tdbs::SOMAGroup::create(sctx, uri, type, {}, tsrng);
}

// [[Rcpp::export]]
Rcpp::XPtr<somaobj_wrap_t> c_group_open(
    std::string& uri,
    std::string& type,
    Rcpp::XPtr<somactx_wrap_t> ctxxp,
    Rcpp::Nullable<Rcpp::DatetimeVector> timestamp = R_NilValue) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // shared pointer to TileDB Context from SOMAContext
    // std::shared_ptr<tiledb::Context> ctx = sctx->tiledb_ctx();

    // optional timestamp range
    std::optional<tdbs::TimestampRange> tsrng = makeTimestampRange(timestamp);

    // Note: both OpenMode.soma_write and OpenMode.soma_delete should be opened
    // in TILEDB_WRITE mode.
    OpenMode mode = type == "READ" ? OpenMode::soma_read : OpenMode::soma_write;

    std::shared_ptr<tdbs::SOMAObject> sgrpptr = tdbs::SOMAObject::open(uri, mode, sctx, tsrng);

    return make_xptr<somaobj_wrap_t>(new SOMAWrapper(sgrpptr));
}

// [[Rcpp::export]]
double c_group_member_count(Rcpp::XPtr<somaobj_wrap_t> xp) {
    check_xptr_tag<somaobj_wrap_t>(xp);  // throws if mismatched
    // unique pointer to SOMAGroup from external pointer wrapper
    return static_cast<double>(xp->ptr<tdbs::SOMAGroup>()->count());
}

// [[Rcpp::export]]
bool c_group_has_member(Rcpp::XPtr<somaobj_wrap_t> xp, const std::string& key) {
    check_xptr_tag<somaobj_wrap_t>(xp);  // throws if mismatched

    return xp->ptr<tdbs::SOMAGroup>()->has(key);
}

// [[Rcpp::export]]
Rcpp::XPtr<somaobj_wrap_t> c_collection_get_member(Rcpp::XPtr<somaobj_wrap_t> xp, const std::string& key) {
    check_xptr_tag<somaobj_wrap_t>(xp);  // throws if mismatched

    somaobj_wrap_t* somaobj_p = new SOMAWrapper(xp->ptr<tdbs::SOMACollectionBase>()->get(key));
    Rcpp::XPtr<somaobj_wrap_t> somaobj_xptr = make_xptr<somaobj_wrap_t>(somaobj_p);

    return somaobj_xptr;
}

// [[Rcpp::export]]
Rcpp::List c_group_members(Rcpp::XPtr<somaobj_wrap_t> xp) {
    check_xptr_tag<somaobj_wrap_t>(xp);  // throws if mismatched
    // unique pointer to SOMAGroup from external pointer wrapper
    // 'members map': a map from string to pair of strings
    auto mm = xp->ptr<tdbs::SOMAGroup>()->members_map();
    auto n = mm.size();
    Rcpp::List lst(n);
    Rcpp::CharacterVector names(n);
    int i = 0;
    for (auto const& key : mm) {
        names[i] = key.first;
        Rcpp::List row = Rcpp::List::create(
            Rcpp::Named("type") = std::string(key.second.second),
            Rcpp::Named("uri") = std::string(key.second.first),
            Rcpp::Named("name") = key.first);
        lst[i] = row;
        i++;
    }
    lst.attr("names") = names;
    return lst;
}

// [[Rcpp::export]]
Rcpp::List c_group_get_metadata(Rcpp::XPtr<somaobj_wrap_t> xp) {
    check_xptr_tag<somaobj_wrap_t>(xp);  // throws if mismatched
    return soma_get_all_metadata(xp);
}

// [[Rcpp::export]]
void c_group_close(Rcpp::XPtr<somaobj_wrap_t> xp) {
    check_xptr_tag<somaobj_wrap_t>(xp);  // throws if mismatched
    // unique pointer to SOMAGroup from external pointer wrapper
    xp->ptr<tdbs::SOMACollectionBase>()->close(true);
}

std::map<int, URIType> uritypemap = {{0, URIType::automatic}, {1, URIType::absolute}, {2, URIType::relative}};

// [[Rcpp::export]]
void c_group_set(
    Rcpp::XPtr<somaobj_wrap_t> xp,
    const std::string& uri,
    int uri_type_int,  // "automatic", "absolute", "relative"
    const std::string& name,
    const std::string& soma_type) {
    check_xptr_tag<somaobj_wrap_t>(xp);  // throws if mismatched
    xp->ptr<tdbs::SOMACollectionBase>()->set(uri, uritypemap[uri_type_int], name, soma_type);
}

// [[Rcpp::export]]
void c_group_remove_member(Rcpp::XPtr<somaobj_wrap_t> xp, const std::string& name) {
    check_xptr_tag<somaobj_wrap_t>(xp);  // throws if mismatched
    xp->ptr<tdbs::SOMACollectionBase>()->del(name);
}

// [[Rcpp::export]]
void c_group_put_metadata(Rcpp::XPtr<somaobj_wrap_t> xp, std::string key, SEXP obj) {
    check_xptr_tag<somaobj_wrap_t>(xp);  // throws if mismatched
    // we implement a simpler interface here as the 'type' is given from the
    // supplied SEXP, as is the extent

    soma_set_metadata(xp, key, obj);
}
