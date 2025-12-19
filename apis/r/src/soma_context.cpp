#include <Rcpp/Lighter>  // for R interface to C++

#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"  // local declarations
#include "xptr-utils.h"  // xptr taggging utilities

// [[Rcpp::export]]
Rcpp::XPtr<somactx_wrap_t> create_soma_context(Rcpp::Nullable<Rcpp::CharacterVector> config = R_NilValue) {
    std::shared_ptr<tiledbsoma::SOMAContext> somactx;
    if (config.isNotNull()) {
        std::map<std::string, std::string> smap;
        auto config_vec = config.as();
        auto config_names = Rcpp::as<Rcpp::CharacterVector>(config_vec.names());
        for (auto& name : config_names) {
            std::string param = Rcpp::as<std::string>(name);
            std::string value = Rcpp::as<std::string>(config_vec[param]);
            smap[param] = value;
        }
        somactx = std::make_shared<tiledbsoma::SOMAContext>(smap);
    } else {
        somactx = std::make_shared<tiledbsoma::SOMAContext>();
    }

    auto ptr = new somactx_wrap_t(somactx);
    auto xp = make_xptr<somactx_wrap_t>(ptr);
    return xp;
}

// [[Rcpp::export]]
Rcpp::CharacterVector get_config_from_soma_context(Rcpp::XPtr<somactx_wrap_t> soma_context) {
    auto tiledb_context = soma_context->ctxptr->tiledb_ctx();
    Rcpp::CharacterVector keys{};
    Rcpp::CharacterVector values{};
    for (const auto& [k, v] : tiledb_context->config()) {
        keys.push_back(k);
        values.push_back(v);
    }
    values.attr("names") = keys;
    return values;
}

// [[Rcpp::export]]
std::string get_data_protocol_from_soma_context(Rcpp::XPtr<somactx_wrap_t> soma_context, const std::string& uri) {
    return soma_context->ctxptr->data_protocol(uri);
}
