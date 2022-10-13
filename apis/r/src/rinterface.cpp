#include <Rcpp.h>
#include <tiledbsoma/tiledbsoma>
#include "pyapi/arrow_adapter.h"
#include <archAPI.h>
#include "rutilities.h"

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

//' @rdname get_table
//' @export
// [[Rcpp::export]]
Rcpp::List export_arrow_array(const std::string& uri,
                              const std::vector<std::string>& colnames,
                              Rcpp::Nullable<Rcpp::XPtr<tiledb::QueryCondition>> qc = R_NilValue,
                              Rcpp::Nullable<Rcpp::DataFrame> dim_points = R_NilValue,
                              Rcpp::Nullable<Rcpp::DataFrame> dim_ranges = R_NilValue,
                              const std::string& loglevel = "warn") {

    tdbs::LOG_SET_LEVEL(loglevel);

    tdbs::LOG_INFO(fmt::format("Reading from {}", uri));

    // Read selected columns from the obs array (obs is unique_ptr<SOMAReader>)
    auto obs = tdbs::SOMAReader::open(uri, "", {}, colnames);

    // If we have a query condition, apply it
    if (!qc.isNull()) {
        tdbs::LOG_INFO(fmt::format("Applying query condition"));
        Rcpp::XPtr<tiledb::QueryCondition> qcxp(qc);
        obs->set_condition(*qcxp);
    }

    // If we have a dimension points, apply them
    if (!dim_points.isNull()) {
        Rcpp::DataFrame df(dim_points);
        Rcpp::CharacterVector nm = df[0];
        Rcpp::CharacterVector tp = df[1];
        Rcpp::NumericVector val = df[2]; // works as proxy for int and float types not for string
        for (int i=0; i<df.nrow(); i++) {
            std::string s(nm[i]);
            std::string t(tp[i]);
            if (t == "UINT64") {
                uint64_t v = static_cast<uint64_t>(makeScalarInteger64(val[i]));
                obs->set_dim_point<uint64_t>(s, v);
                tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} for {} with {}", i, nm[i], tp[i], v));
            } else if (t == "INT64") {
                int64_t v = makeScalarInteger64(val[i]);
                obs->set_dim_point<int64_t>(s, v);
                tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} for {} with {}", i, nm[i], tp[i], v));
            } else if (t == "FLOAT32") {
                float v = static_cast<float>(val[i]);
                obs->set_dim_point<float>(s, v);
                tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} for {} with {}", i, nm[i], tp[i], v));
            } else if (t == "FLOAT64") {
                double v = val[i];
                obs->set_dim_point<double>(s, v);
                tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} for {} with {}", i, nm[i], tp[i], v));
            } else if (t == "INT32") {
                int32_t v = static_cast<int32_t>(val[i]);
                obs->set_dim_point<int32_t>(s, v);
                tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} for {} with {}", i, nm[i], tp[i], v));
            } else {
                Rcpp::stop("Currently unsupported type: ", t);
            }
        }
    }

    // If we have a dimension points, apply them
    if (!dim_ranges.isNull()) {
        Rcpp::DataFrame df(dim_ranges);
        Rcpp::CharacterVector nm = df[0];
        Rcpp::CharacterVector tp = df[1];
        Rcpp::NumericVector lo = df[2]; // works as proxy for int and float types not for string
        Rcpp::NumericVector hi = df[3]; // works as proxy for int and float types not for string
        for (int i=0; i<df.nrow(); i++) {
            std::string s(nm[i]);
            std::string t(tp[i]);
            if (t == "UINT64") {
                uint64_t l = static_cast<uint64_t>(makeScalarInteger64(lo[i]));
                uint64_t h = static_cast<uint64_t>(makeScalarInteger64(hi[i]));
                std::vector<std::pair<uint64_t, uint64_t>> vp{std::make_pair(l,h)};
                obs->set_dim_ranges<uint64_t>(s, vp);
                tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} for {} with {} - {}", i, nm[i], tp[i], l, h));
            } else if (t == "INT64") {
                int64_t l = makeScalarInteger64(lo[i]);
                int64_t h = makeScalarInteger64(hi[i]);
                std::vector<std::pair<int64_t, int64_t>> vp{std::make_pair(l,h)};
                obs->set_dim_ranges<int64_t>(s, vp);
                tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} for {} with {} - {}", i, nm[i], tp[i], l, h));
            } else if (t == "FLOAT32") {
                float l = static_cast<float>(lo[i]);
                float h = static_cast<float>(hi[i]);
                std::vector<std::pair<float, float>> vp{std::make_pair(l,h)};
                obs->set_dim_ranges<float>(s, vp);
                tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} for {} with {} - {}", i, nm[i], tp[i], l, h));
            } else if (t == "FLOAT64") {
                double l = lo[i];
                double h = hi[i];
                std::vector<std::pair<double, double>> vp{std::make_pair(l,h)};
                obs->set_dim_ranges<double>(s, vp);
                tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} for {} with {} - {}", i, nm[i], tp[i], l, h));
            } else if (t == "INT32") {
                int32_t l = static_cast<int32_t>(lo[i]);
                int32_t h = static_cast<int32_t>(hi[i]);
                std::vector<std::pair<int32_t, int32_t>> vp{std::make_pair(l,h)};
                obs->set_dim_ranges<int32_t>(s, vp);
                tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} for {} with {} - {}", i, nm[i], tp[i], l, h));
            } else {
                Rcpp::stop("Currently unsupported type: ", t);
            }
        }
    }

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
    Rcpp::List as = Rcpp::List::create(Rcpp::Named("schema") = sxp,
                                       Rcpp::Named("array_data") = axp);
    as.attr("class") = "arch_array";
    return as;
}

//' @noRd
// [[Rcpp::export]]
void set_log_level(const std::string& level) {
    tdbs::LOG_SET_LEVEL(level);
}

//' @rdname get_table
//' @export
// [[Rcpp::export]]
double nnz(const std::string& uri) {
    auto sr = tdbs::SOMAReader::open(uri);
    return static_cast<double>(sr->nnz());
}

//' @rdname get_table
//' @export
// [[Rcpp::export]]
Rcpp::CharacterVector get_column_types(const std::string& uri,
                                       const std::vector<std::string>& names) {

    auto obs = tdbs::SOMAReader::open(uri);
    obs->submit();
    auto obs_data = obs->read_next();
    size_t n = names.size();
    Rcpp::CharacterVector vs(n);
    for (size_t i=0; i<n; i++) {
        auto datatype = obs_data->get()->at(names[i])->type();
        vs[i] = std::string(_tiledb_datatype_to_string(datatype));
    }
    vs.attr("names") = names;
    return vs;
}
