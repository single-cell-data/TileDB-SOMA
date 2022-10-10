#include <Rcpp.h>
#include <tiledbsoma/tiledbsoma>
#include "pyapi/arrow_adapter.h"
#include <archAPI.h>

namespace tdbs = tiledbsoma;

// initial 'proof of light' function to access column names
// will likely be replaced / changed in due course
//
//' @rdname get_table
//' @export
// [[Rcpp::export]]
std::vector<std::string> get_column_names(const std::string& uri) {
    std::map<std::string, std::string> config;

    // Read all values from the obs array (obs is unique_ptr<SOMAReader>)
    auto obs = tdbs::SOMAReader::open(uri);
    obs->submit();

    // Getting next batch:  std::optinal<std::shared_ptr<ArrayBuffers>>
    auto obs_data = obs->read_next();
    if (!obs->results_complete()) {
        Rcpp::warning("Read of '%s' incomplete", uri);
    }

    // names gets a vector of strings, at() gets a particular column buffer by name
    auto names = obs_data->get()->names();

    return names;
}

// initial 'proof of light' function to a column by name
// will likely be replaced / changed in due course
//
// [[Rcpp::export]]
bool export_column(const std::string& uri, const std::string& colname,
                   SEXP schemaxp, SEXP arrayxp) {
    std::map<std::string, std::string> config;

    // Read all values from the obs array (obs is unique_ptr<SOMAReader>)
    auto obs = tdbs::SOMAReader::open(uri);
    obs->submit();

    // Getting next batch:  std::optinal<std::shared_ptr<ArrayBuffers>>
    auto obs_data = obs->read_next();
    if (!obs->results_complete()) {
        Rcpp::warning("Read of '%s' incomplete", uri);
    }

    // now buf is a shared_ptr to ColumnBuffer
    auto buf = obs_data->get()->at(colname);

    // this is pair of array and schema pointer
    auto pp = tdbs::ArrowAdapter::to_arrow(buf);

    memcpy((void*) R_ExternalPtrAddr(schemaxp), pp.second.get(), sizeof(ArrowSchema));
    memcpy((void*) R_ExternalPtrAddr(arrayxp), pp.first.get(), sizeof(ArrowArray));

    return true;
}


// more efficient column accessor doing more 'in C++' rather than
// still relies on helper from package 'arch' that we expect to be
// replaced in due course by package 'nanoarrow' (once released)
//
//' @rdname get_table
//' @export
// [[Rcpp::export]]
SEXP export_column_direct(const std::string& uri, const std::vector<std::string>& colnames) {

    tdbs::LOG_INFO(fmt::format("Reading from {}", uri));

    // Read selected columns from the obs array (obs is unique_ptr<SOMAReader>)
    auto obs = tdbs::SOMAReader::open(uri, "", {}, colnames);
    obs->submit();

    // Getting next batch:  std::optinal<std::shared_ptr<ArrayBuffers>>
    auto obs_data = obs->read_next();
    if (!obs->results_complete()) {
        Rcpp::warning("Read of '%s' incomplete", uri);
    }
    tdbs::LOG_INFO(fmt::format("Read complete with {} obs and {} cols",
                               obs_data->get()->num_rows(),
                               obs_data->get()->names().size()));

    auto ncol = obs_data->get()->names().size();
    Rcpp::List reslist(ncol);
    reslist.attr("names") = obs_data->get()->names();

    for (size_t i=0; i<ncol; i++) {
        // this allocates, and properly wraps as external pointers controlling lifetime
        SEXP schemaxp = arch_c_allocate_schema();
        SEXP arrayxp = arch_c_allocate_array_data();

        // now buf is a shared_ptr to ColumnBuffer
        auto buf = obs_data->get()->at(colnames[i]);

        tdbs::LOG_INFO(fmt::format("Accessing {} at {}", colnames[i], i));

        // this is pair of array and schema pointer
        auto pp = tdbs::ArrowAdapter::to_arrow(buf);

        memcpy((void*) R_ExternalPtrAddr(schemaxp), pp.second.get(), sizeof(ArrowSchema));
        memcpy((void*) R_ExternalPtrAddr(arrayxp), pp.first.get(), sizeof(ArrowArray));

        // in R:
        //   aa <- arch::arch_array(schema, array, FALSE)    # R/array.R:16
        // which just passes throught and then creates as S3 class
        // sp here we just set an S3 object up: a list with a class attribute
        Rcpp::List res = Rcpp::List::create(Rcpp::Named("schema") = schemaxp,
                                            Rcpp::Named("array_data") = arrayxp);
        res.attr("class") = Rcpp::CharacterVector::create("arch_array");
        // in R:  arch::from_arch_array()

        reslist[i] = res;
    }

    return reslist;
}

inline SEXP array_sexp_new(SEXP schema_xptr, SEXP array_data_xptr) {
    const char* names[] = {"schema", "array_data", ""};
    SEXP array_sexp = PROTECT(Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(array_sexp, 0, schema_xptr);
    SET_VECTOR_ELT(array_sexp, 1, array_data_xptr);
    Rf_setAttrib(array_sexp, R_ClassSymbol, Rf_mkString("arch_array"));
    UNPROTECT(1);
    return array_sexp;
}

//' @rdname get_table
//' @export
// [[Rcpp::export]]
SEXP export_recordbatch(const std::string& uri, const std::vector<std::string>& colnames) {

    tdbs::LOG_INFO(fmt::format("Reading from {}", uri));

    // Read selected columns from the obs array (obs is unique_ptr<SOMAReader>)
    auto obs = tdbs::SOMAReader::open(uri, "", {}, colnames);
    obs->submit();

    // Getting next batch:  std::optional<std::shared_ptr<ArrayBuffers>>
    auto obs_data = obs->read_next();
    if (!obs->results_complete()) {
        Rcpp::warning("Read of '%s' incomplete", uri);
    }
    tdbs::LOG_INFO(fmt::format("Read complete with {} obs and {} cols",
                               obs_data->get()->num_rows(),
                               obs_data->get()->names().size()));

    auto ncol = obs_data->get()->names().size();
    Rcpp::List schlst(ncol), arrlst(ncol);

    for (size_t i=0; i<ncol; i++) {
        // this allocates, and properly wraps as external pointers controlling lifetime
        SEXP schemaxp = arch_c_allocate_schema();
        SEXP arrayxp = arch_c_allocate_array_data();

        // now buf is a shared_ptr to ColumnBuffer
        auto buf = obs_data->get()->at(colnames[i]);

        tdbs::LOG_INFO(fmt::format("Accessing {} at {}", colnames[i], i));

        // this is pair of array and schema pointer
        auto pp = tdbs::ArrowAdapter::to_arrow(buf);

        pp.second->name = colnames[i].c_str();

        memcpy((void*) R_ExternalPtrAddr(schemaxp), pp.second.get(), sizeof(ArrowSchema));
        memcpy((void*) R_ExternalPtrAddr(arrayxp), pp.first.get(), sizeof(ArrowArray));

        schlst[i] = schemaxp;
        arrlst[i] = arrayxp;
    }

    struct ArrowArray* array_data_tmp = (struct ArrowArray*) R_ExternalPtrAddr(arrlst[0]);
    int rows = static_cast<int>(array_data_tmp->length);
    SEXP sxp = arch_c_schema_xptr_new(Rcpp::wrap("+s"), 	// format
                                      Rcpp::wrap(""),   	// name
                                      Rcpp::List(),       	// metadata
                                      Rcpp::wrap(2),      	// flags, 2 == unordered, nullable, no sorted map keys
                                      schlst, 	        	// children
                                      R_NilValue);        	// dictionary
    SEXP axp = arch_c_array_from_sexp(Rcpp::List::create(Rcpp::Named("")=R_NilValue), // buffers
                                      Rcpp::wrap(rows), 	// length
                                      Rcpp::wrap(-1), 	    // null count, -1 means not determined
                                      Rcpp::wrap(0),    	// offset (in bytes)
                                      arrlst,               // children
                                      R_NilValue);          // dictionary
    SEXP as = array_sexp_new(sxp, axp);
    return as;
}
