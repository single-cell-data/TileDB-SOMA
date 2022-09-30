
#include <Rcpp.h>
#include <tiledbsoma/tiledbsoma>
#include "pyapi/arrow_adapter.h"

namespace tdbs = tiledbsoma;

// Illustrative first example, converted from cli.cc
// // [[Rcpp::export]]
// bool read_sdf(const std::string& uri) {
//     std::map<std::string, std::string> config;

//     // Read all values from the obs array
//     // obs is unique_ptr<SOMAReader>
//     auto obs = tdbs::SOMAReader::open(uri + "/obs", "obs");
//     obs->submit();
//     // Getting next batch:  std::optinal<std::shared_ptr<ArrayBuffers>>
//     auto obs_data = obs->read_next();

//     // names gets a vector of strings, at() gets a particular column buffer by name
//     auto names = obs_data->get()->names();
//     for (int i=0; i<names.size(); i++) {
//         // now buf is a shared_ptr to ColumnBuffer
//         auto buf = obs_data->get()->at(names[i]);
//         Rcpp::Rcout << names[i] << " " << buf->type() << std::endl;

//         // this is pair of array and schema pointer
//         //auto pp = tdbs::ArrowAdapter::to_arrow(buf);

//     }
//     return obs->results_complete();
// }

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

// // [[Rcpp::export]]
// Rcpp::List get_column(const std::string& uri, const std::string& name) {
//     std::map<std::string, std::string> config;

//     // Read all values from the obs array
//     // obs is unique_ptr<SOMAReader>
//     auto obs = tdbs::SOMAReader::open(uri + "/obs", "obs");
//     obs->submit();
//     // Getting next batch:  std::optinal<std::shared_ptr<ArrayBuffers>>
//     auto obs_data = obs->read_next();

//     // now buf is a shared_ptr to ColumnBuffer
//     auto buf = obs_data->get()->at(name);

//     // this is pair of array and schema pointer
//     auto pp = tdbs::ArrowAdapter::to_arrow(buf);

//     Rcpp::XPtr<ArrowArray> ap = Rcpp::XPtr<ArrowArray>(pp.first.get());
//     Rcpp::XPtr<ArrowSchema> sp = Rcpp::XPtr<ArrowSchema>(pp.second.get());
//     ArrowArray a = *ap;
//     ArrowSchema s = *sp;
//     Rcpp::Rcout << "length: " << a.length << std::endl
//                 << "nbuff: " << a.n_buffers << std::endl
//                 << "format: " << s.format << std::endl;

//     //ArrowArray a = *pp.first.get();
//     //ArrowSchema s = *pp.second.get();

//     return Rcpp::List::create(Rcpp::Named("array") = ap,
//                               Rcpp::Named("schema") = sp);
// }

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
