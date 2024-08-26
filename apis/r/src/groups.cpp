#include <Rcpp/Lighter>                         // for R interface to C++
#include <nanoarrow/r.h>                        // for C/C++ interface to Arrow (via header exported from the R package)
#include <nanoarrow/nanoarrow.hpp>              // for C/C++ interface to Arrow (vendored)
#include <RcppInt64>                            // for fromInteger64

#include <tiledbsoma/tiledbsoma>
#include <tiledbsoma/reindexer/reindexer.h>

#include "rutilities.h"         				// local declarations
#include "xptr-utils.h"         				// xptr taggging utilities

namespace tdbs = tiledbsoma;

//' Create a group
//' @param uri The array URI
//' @param ctxxp An external pointer to the SOMAContext wrapper
//' @export
// [[Rcpp::export]]
void c_group_create(std::string& uri, std::string& type, Rcpp::XPtr<somactx_wrap_t> ctxxp,
                    Rcpp::Nullable<Rcpp::DatetimeVector> timestamp = R_NilValue) {

    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // shared pointer to TileDB Context from SOMAContext
    // not needed here:  std::shared_ptr<tiledb::Context> ctx = sctx->tiledb_ctx();

    // optional timestamp range
    std::optional<tdbs::TimestampRange> tsrng = makeTimestampRange(timestamp);
    if (timestamp.isNotNull()) {
        Rcpp::DatetimeVector v(timestamp);
        spdl::debug("[c_group_create] uri {} ts ({},{})", uri, v[0], v[1]);
    } else {
        spdl::debug("[c_group_create] uri {}", uri);
    }

    tdbs::SOMAGroup::create(sctx, uri, type, tsrng);
}

// [[Rcpp::export]]
Rcpp::XPtr<somagrp_wrap_t> c_group_open(std::string& uri, std::string& type,
                                        Rcpp::XPtr<somactx_wrap_t> ctxxp,
                                        Rcpp::Nullable<Rcpp::DatetimeVector> timestamp = R_NilValue) {

    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // shared pointer to TileDB Context from SOMAContext
    // std::shared_ptr<tiledb::Context> ctx = sctx->tiledb_ctx();

    // optional timestamp range
    std::optional<tdbs::TimestampRange> tsrng = makeTimestampRange(timestamp);

    OpenMode mode = type == "READ" ? OpenMode::read : OpenMode::write;

    auto sgrpptr = tdbs::SOMAGroup::open(mode, uri, sctx, "unnamed", tsrng);

    somagrp_wrap_t* somagrp_p = new SOMAGroupWrapper(std::move(sgrpptr));
    Rcpp::XPtr<somagrp_wrap_t> somagrp_xptr = make_xptr<somagrp_wrap_t>(somagrp_p);
    return somagrp_xptr;
}


// [[Rcpp::export]]
double c_group_member_count(Rcpp::XPtr<somagrp_wrap_t> xp) {
    check_xptr_tag<somagrp_wrap_t>(xp);  // throws if mismatched
    // unique pointer to SOMAGroup from external pointer wrapper
    return static_cast<double>( xp->grpptr->count() );
}

// [[Rcpp::export]]
Rcpp::List c_group_members(Rcpp::XPtr<somagrp_wrap_t> xp) {
    check_xptr_tag<somagrp_wrap_t>(xp);  // throws if mismatched
    // unique pointer to SOMAGroup from external pointer wrapper
    // 'members map': a map from string to pair of strings
    auto mm = xp->grpptr->members_map();
    auto n = mm.size();
    Rcpp::List lst(n);
    Rcpp::CharacterVector names(n);
    int i = 0;
    for (auto const& key: mm) {
        names[i] = key.first;
        Rcpp::List row = Rcpp::List::create(Rcpp::Named("type") = std::string(key.second.second),
                                            Rcpp::Named("uri") = std::string(key.second.first),
                                            Rcpp::Named("name") = key.first);
        lst[i] = row;
        i++;
    }
    lst.attr("names") = names;
    return lst;
}

// borrowed with a tip-of-the-hat from tiledb::src/libtiledb.coo
// helper function to copy int vector
template <typename T>
Rcpp::IntegerVector copy_int_vector(const uint32_t v_num, const void* v) {
  // Strictly speaking a check for under/overflow would be needed here yet this for
  // metadata annotation (and not data payload) so extreme ranges are less likely
  Rcpp::IntegerVector vec(v_num);
  const T *ivec = static_cast<const T*>(v);
  size_t n = static_cast<size_t>(v_num);
  for (size_t i=0; i<n; i++) vec[i] = static_cast<int32_t>(ivec[i]);
  return(vec);
}
// helper function to convert_metadata
SEXP _metadata_to_sexp(const tiledb_datatype_t v_type, const uint32_t v_num, const void* v) {
    // This supports a limited set of basic types as the metadata
    // annotation is not meant to support complete serialization
    if (v_type == TILEDB_INT32) {
        Rcpp::IntegerVector vec(v_num);
        std::memcpy(vec.begin(), v, v_num*sizeof(int32_t));
        return(vec);
    } else if (v_type == TILEDB_FLOAT64) {
        Rcpp::NumericVector vec(v_num);
        std::memcpy(vec.begin(), v, v_num*sizeof(double));
        return(vec);
    } else if (v_type == TILEDB_FLOAT32) {
        Rcpp::NumericVector vec(v_num);
        const float *fvec = static_cast<const float*>(v);
        size_t n = static_cast<size_t>(v_num);
        for (size_t i=0; i<n; i++) vec[i] = static_cast<double>(fvec[i]);
        return(vec);
    } else if (v_type == TILEDB_CHAR || v_type == TILEDB_STRING_ASCII || v_type == TILEDB_STRING_UTF8) {
        std::string s(static_cast<const char*>(v), v_num);
        return(Rcpp::wrap(s));
    } else if (v_type == TILEDB_INT8) {
        Rcpp::LogicalVector vec(v_num);
        const int8_t *ivec = static_cast<const int8_t*>(v);
        size_t n = static_cast<size_t>(v_num);
        for (size_t i=0; i<n; i++) vec[i] = static_cast<bool>(ivec[i]);
        return(vec);
    } else if (v_type == TILEDB_UINT8) {
        // Strictly speaking a check for under/overflow would be needed here (and below) yet this
        // is for metadata annotation (and not data payload) so extreme ranges are less likely
        return copy_int_vector<uint8_t>(v_num, v);
    } else if (v_type == TILEDB_INT16) {
        return copy_int_vector<int16_t>(v_num, v);
    } else if (v_type == TILEDB_UINT16) {
        return copy_int_vector<uint16_t>(v_num, v);
    } else if (v_type == TILEDB_UINT32) {
        return copy_int_vector<uint32_t>(v_num, v);
    } else if (v_type == TILEDB_INT64) {
        std::vector<int64_t> iv(v_num);
        std::memcpy(&(iv[0]), v, v_num*sizeof(int64_t));
        return Rcpp::toInteger64(iv);
    } else if (v_type == TILEDB_UINT64) {
        return copy_int_vector<uint64_t>(v_num, v);
    } else {
        Rcpp::stop("No support yet for TileDB data type %s", tiledb::impl::type_to_str(v_type));
  }
}



// [[Rcpp::export]]
Rcpp::List c_group_get_metadata(Rcpp::XPtr<somagrp_wrap_t> xp) {
    check_xptr_tag<somagrp_wrap_t>(xp);  // throws if mismatched
    // unique pointer to SOMAGroup from external pointer wrapper
    // 'members map': a map from string to pair of strings
    auto mm = xp->grpptr->get_metadata();
    auto n = mm.size();
    Rcpp::List lst(n);
    Rcpp::CharacterVector names(n);
    int i = 0;
    for (auto const& key: mm) {
        names[i] = key.first;
        // what is keyed using MetadataValue = std::tuple<tiledb_datatype_t, uint32_t, const void*>
        auto tpl = key.second;
        auto sxp = _metadata_to_sexp(std::get<0>(tpl), std::get<1>(tpl), std::get<2>(tpl));
        Rcpp::List row = Rcpp::List::create(Rcpp::Named("name") = sxp);
        lst[i] = row;
        i++;
    }
    lst.attr("names") = names;
    return lst;
}

// [[Rcpp::export]]
void c_group_close(Rcpp::XPtr<somagrp_wrap_t> xp) {
    check_xptr_tag<somagrp_wrap_t>(xp);  // throws if mismatched
    // unique pointer to SOMAGroup from external pointer wrapper
    xp->grpptr->close();
}

std::map<int, URIType> uritypemap = { { 0, URIType::automatic },
                                      { 1, URIType::absolute },
                                      { 2, URIType::relative } };


// [[Rcpp::export]]
void c_group_set(Rcpp::XPtr<somagrp_wrap_t> xp,
                 const std::string& uri,
                 int uri_type_int,  							// "automatic", "absolute", "relative"
                 const std::string& name,
                 const std::string& soma_type) {
    check_xptr_tag<somagrp_wrap_t>(xp);  // throws if mismatched
    xp->grpptr->set(uri, uritypemap[uri_type_int], name, soma_type);
}

// [[Rcpp::export]]
void c_group_remove_member(Rcpp::XPtr<somagrp_wrap_t> xp, const std::string& name) {
    check_xptr_tag<somagrp_wrap_t>(xp);  // throws if mismatched
    xp->grpptr->del(name);
}

// [[Rcpp::export]]
void c_group_put_metadata(Rcpp::XPtr<somagrp_wrap_t> xp, std::string key, SEXP obj) {
    check_xptr_tag<somagrp_wrap_t>(xp);  // throws if mismatched
    // we implement a simpler interface here as the 'type' is given from the supplied SEXP, as is the extent
    switch(TYPEOF(obj)) {
    case VECSXP: {
        Rcpp::stop("List objects are not supported.");
        break;// not reached
        }
    case REALSXP: {
        Rcpp::NumericVector v(obj);
        if (Rcpp::isInteger64(obj)) {
            xp->grpptr->set_metadata(key, TILEDB_INT64, v.size(), v.begin());
        } else {
            xp->grpptr->set_metadata(key, TILEDB_FLOAT64, v.size(), v.begin());
        }
        break;
    }
    case INTSXP: {
        Rcpp::IntegerVector v(obj);
        xp->grpptr->set_metadata(key, TILEDB_INT32, v.size(), v.begin());
        break;
    }
    case STRSXP: {
        Rcpp::CharacterVector v(obj);
        std::string s(v[0]);
        // We use TILEDB_CHAR interchangeably with TILEDB_STRING_ASCII is this best string type?
        xp->grpptr->set_metadata(key, TILEDB_STRING_ASCII, s.length(), s.c_str());
        break;
    }
    case LGLSXP: {              // experimental: map R logical (ie TRUE, FALSE, NA) to int8
        Rcpp::stop("Logical vector objects are not supported.");
        break;// not reached
    }
    default: {
        Rcpp::stop("No support (yet) for type '%s'.", Rf_type2char(TYPEOF(obj)));
        break; // not reached
    }
    }

}
