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
                              Rcpp::Nullable<Rcpp::List> dim_points = R_NilValue,
                              Rcpp::Nullable<Rcpp::List> dim_ranges = R_NilValue,
                              const std::string& loglevel = "warn") {

    tdbs::LOG_SET_LEVEL(loglevel);

    tdbs::LOG_INFO(fmt::format("Reading from {}", uri));

    // Read selected columns from the obs array (obs is unique_ptr<SOMAReader>)
    auto obs = tdbs::SOMAReader::open(uri, "", {}, colnames);

    std::unordered_map<std::string, tiledb_datatype_t> name2type;

    std::shared_ptr<tiledb::ArraySchema> schema = obs->schema();
    tiledb::Domain domain = schema->domain();
    std::vector<tiledb::Dimension> dims = domain.dimensions();
    for (auto& dim: dims) {
        tdbs::LOG_INFO(fmt::format("Dimension {} type {} domain {} extent {}",
                                   dim.name(), tiledb::impl::to_str(dim.type()),
                                   dim.domain_to_str(), dim.tile_extent_to_str()));
        name2type.emplace(std::make_pair(dim.name(), dim.type()));
    }

    // If we have a query condition, apply it
    if (!qc.isNull()) {
        tdbs::LOG_INFO(fmt::format("Applying query condition"));
        Rcpp::XPtr<tiledb::QueryCondition> qcxp(qc);
        obs->set_condition(*qcxp);
    }

    // If we have dimension points, apply them
    // The interface is named list, where each (named) list elements is one (named) dimension
    // The List element is a simple vector of points and each point is applied to the named dimension
    if (!dim_points.isNull()) {
        Rcpp::List lst(dim_points);
        std::vector<std::string> colnames = lst.attr("names");
        for (auto& nm: colnames) {
            auto tp = name2type[nm];
            if (tp == TILEDB_UINT64) {
                Rcpp::NumericVector payload = lst[nm];
                std::vector<int64_t> iv = getInt64Vector(payload);
                std::vector<uint64_t> uv(iv.size());
                for (size_t i=0; i<iv.size(); i++) {
                    uv[i] = static_cast<uint64_t>(iv[i]);
                    obs->set_dim_point<uint64_t>(nm, uv[i]);  // bonked when use with vector
                    tdbs::LOG_INFO(fmt::format("Applying dim point {} on {}", uv[i], nm));
                }
            } else if (tp == TILEDB_INT64) {
                Rcpp::NumericVector payload = lst[nm];
                std::vector<int64_t> iv = getInt64Vector(payload);
                for (size_t i=0; i<iv.size(); i++) {
                    obs->set_dim_point<int64_t>(nm, iv[i]);
                    tdbs::LOG_INFO(fmt::format("Applying dim point {} on {}", iv[i], nm));
                }
            } else if (tp == TILEDB_FLOAT32) {
                Rcpp::NumericVector payload = lst[nm];
                for (R_xlen_t i=0; i<payload.size(); i++) {
                    float v = static_cast<uint64_t>(payload[i]);
                    obs->set_dim_point<float>(nm, v);
                    tdbs::LOG_INFO(fmt::format("Applying dim point {} on {}", v, nm));
                }
            } else if (tp == TILEDB_FLOAT64) {
                Rcpp::NumericVector payload = lst[nm];
                for (R_xlen_t i=0; i<payload.size(); i++) {
                    obs->set_dim_point<double>(nm,payload[i]);
                    tdbs::LOG_INFO(fmt::format("Applying dim point {} on {}", payload[i], nm));
                }
            } else if (tp == TILEDB_INT32) {
                Rcpp::IntegerVector payload = lst[nm];
                for (R_xlen_t i=0; i<payload.size(); i++) {
                    obs->set_dim_point<int32_t>(nm,payload[i]);
                    tdbs::LOG_INFO(fmt::format("Applying dim point {} on {}", payload[i], nm));
                }
            } else {
                Rcpp::stop("Currently unsupported type: ", tiledb::impl::to_str(tp));
            }
        }
    }

    // If we have a dimension points, apply them
    if (!dim_ranges.isNull()) {
        Rcpp::List lst(dim_ranges);
        std::vector<std::string> colnames = lst.attr("names");
        for (auto& nm: colnames) {
            auto tp = name2type[nm];
            if (tp == TILEDB_UINT64) {
                Rcpp::NumericMatrix mm = lst[nm];
                Rcpp::NumericMatrix::Column lo = mm.column(0); // works as proxy for int and float types
                Rcpp::NumericMatrix::Column hi = mm.column(1); // works as proxy for int and float types
                for (int i=0; i<mm.nrow(); i++) {
                    uint64_t l = static_cast<uint64_t>(makeScalarInteger64(lo[i]));
                    uint64_t h = static_cast<uint64_t>(makeScalarInteger64(hi[i]));
                    std::vector<std::pair<uint64_t, uint64_t>> vp{std::make_pair(l,h)};
                    obs->set_dim_ranges<uint64_t>(nm, vp);
                    tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} with {} - {}", i, nm, l, h));
                }
            } else if (tp == TILEDB_INT64) {
                Rcpp::NumericMatrix mm = lst[nm];
                Rcpp::NumericMatrix::Column lo = mm.column(0); // works as proxy for int and float types
                Rcpp::NumericMatrix::Column hi = mm.column(1); // works as proxy for int and float types
                for (int i=0; i<mm.nrow(); i++) {
                    std::vector<std::pair<int64_t, int64_t>> vp{std::make_pair(lo[i], hi[i])};
                    obs->set_dim_ranges<int64_t>(nm, vp);
                    tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} with {} - {}", i, nm, lo[i], hi[i]));
                }
            } else if (tp == TILEDB_FLOAT32) {
                Rcpp::NumericMatrix mm = lst[nm];
                Rcpp::NumericMatrix::Column lo = mm.column(0); // works as proxy for int and float types
                Rcpp::NumericMatrix::Column hi = mm.column(1); // works as proxy for int and float types
                for (int i=0; i<mm.nrow(); i++) {
                    float l = static_cast<float_t>(lo[i]);
                    float h = static_cast<float_t>(hi[i]);
                    std::vector<std::pair<float, float>> vp{std::make_pair(l,h)};
                    obs->set_dim_ranges<float>(nm, vp);
                    tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} with {} - {}", i, nm, l, h));
                }
            } else if (tp == TILEDB_FLOAT64) {
                Rcpp::NumericMatrix mm = lst[nm];
                Rcpp::NumericMatrix::Column lo = mm.column(0); // works as proxy for int and float types
                Rcpp::NumericMatrix::Column hi = mm.column(1); // works as proxy for int and float types
                for (int i=0; i<mm.nrow(); i++) {
                    std::vector<std::pair<double, double>> vp{std::make_pair(lo[i],hi[i])};
                    obs->set_dim_ranges<double>(nm, vp);
                    tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} with {} - {}", i, nm, lo[i], hi[i]));
                }
            } else if (tp == TILEDB_INT32) {
                Rcpp::IntegerMatrix mm = lst[nm];
                Rcpp::IntegerMatrix::Column lo = mm.column(0); // works as proxy for int and float types
                Rcpp::IntegerMatrix::Column hi = mm.column(1); // works as proxy for int and float types
                for (int i=0; i<mm.nrow(); i++) {
                    int32_t l = static_cast<int32_t>(lo[i]);
                    int32_t h = static_cast<int32_t>(hi[i]);
                    std::vector<std::pair<int32_t, int32_t>> vp{std::make_pair(l,h)};
                    obs->set_dim_ranges<int32_t>(nm, vp);
                    tdbs::LOG_INFO(fmt::format("Applying dim point {} on {} with {} - {}", i, nm[i], l, h));
                }
            } else {
                Rcpp::stop("Currently unsupported type: ", tiledb::impl::to_str(tp));
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
                                       const std::vector<std::string>& colnames) {

    auto obs = tdbs::SOMAReader::open(uri);
    obs->submit();
    auto obs_data = obs->read_next();
    size_t n = colnames.size();
    Rcpp::CharacterVector vs(n);
    for (size_t i=0; i<n; i++) {
        auto datatype = obs_data->get()->at(colnames[i])->type();
        vs[i] = std::string(tiledb::impl::to_str(datatype));
    }
    vs.attr("names") = colnames;
    return vs;
}
