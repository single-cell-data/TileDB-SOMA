#include <Rcpp/Lighter>  // for R interface to C++
#include <sstream>

#include <nanoarrow/r.h>            // for C/C++ interface to Arrow (via header exported from the R package)
#include <RcppInt64>                // for fromInteger64
#include <nanoarrow/nanoarrow.hpp>  // for C/C++ interface to Arrow (vendored)

#include <tiledbsoma/reindexer/reindexer.h>
#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"  // local declarations
#include "xptr-utils.h"  // xptr taggging utilities

namespace tdbs = tiledbsoma;

// [[Rcpp::export]]
Rcpp::List soma_group_get_members(Rcpp::XPtr<tiledbsoma::SOMAGroup> group) {
    auto mm = group->members_map();
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

std::map<int, URIType> uritypemap = {{0, URIType::automatic}, {1, URIType::absolute}, {2, URIType::relative}};

// [[Rcpp::export]]
void soma_group_set(
    Rcpp::XPtr<tiledbsoma::SOMAGroup> group,
    const std::string& uri,
    int uri_type_int,  // "automatic", "absolute", "relative"
    const std::string& name,
    const std::string& soma_type) {
    if (!group) {
        throw Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    group->set(uri, uritypemap[uri_type_int], name, soma_type);
}

// [[Rcpp::export]]
void soma_group_remove_member(Rcpp::XPtr<tiledbsoma::SOMAGroup> group, const std::string& name) {
    if (!group) {
        throw Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    group->del(name);
}
