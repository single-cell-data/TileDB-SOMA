#include <Rcpp.h>
#include <tiledbsoma/tiledbsoma>
#include "pyapi/arrow_adapter.h"

namespace tdbs = tiledbsoma;

// [[Rcpp::export]]
std::vector<std::string> get_column_names(const std::string& uri) {
    std::map<std::string, std::string> config;

    // Read all values from the obs array
    // obs is unique_ptr<SOMAReader>
    auto obs = tdbs::SOMAReader::open(uri);
    obs->submit();
    // Getting next batch:  std::optinal<std::shared_ptr<ArrayBuffers>>
    auto obs_data = obs->read_next();

    // names gets a vector of strings, at() gets a particular column buffer by name
    auto names = obs_data->get()->names();

    return names;
}

// [[Rcpp::export]]
bool export_column(const std::string& uri, const std::string& name, SEXP schemaxp, SEXP arrayxp) {
    std::map<std::string, std::string> config;

    // Read all values from the obs array
    // obs is unique_ptr<SOMAReader>
    auto obs = tdbs::SOMAReader::open(uri);
    obs->submit();
    // Getting next batch:  std::optinal<std::shared_ptr<ArrayBuffers>>
    auto obs_data = obs->read_next();

    // now buf is a shared_ptr to ColumnBuffer
    auto buf = obs_data->get()->at(name);

    // this is pair of array and schema pointer
    auto pp = tdbs::ArrowAdapter::to_arrow(buf);

    memcpy((void*) R_ExternalPtrAddr(schemaxp), pp.second.get(), sizeof(ArrowSchema));
    memcpy((void*) R_ExternalPtrAddr(arrayxp), pp.first.get(), sizeof(ArrowArray));

    return true;
}
