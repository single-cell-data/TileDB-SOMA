#include <Rcpp/Lighter>  // for R interface to C++

#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"  // local declarations
#include "xptr-utils.h"  // xptr taggging utilities

// [[Rcpp::export]]
Rcpp::XPtr<somactx_wrap_t> create_soma_context(Rcpp::Nullable<Rcpp::CharacterVector> config = R_NilValue) {
    // if we hae a config, use it
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
