#include <Rcpp.h>
#include <tiledbsoma/tiledbsoma>
#include <archAPI.h>
#include "rutilities.h"

namespace tdbs = tiledbsoma;

// Initial 'proof of light' function to access column names
// will likely be replaced / changed in due course
//
//' @noRd
// [[Rcpp::export]]
std::vector<std::string> get_column_names(const std::string& uri) {
    std::map<std::string, std::string> config;

    // Read all values from the array (return is unique_ptr<SOMAReader>)
    auto sr = tdbs::SOMAReader::open(uri);
    sr->submit();

    // Getting next batch:  std::optinal<std::shared_ptr<ArrayBuffers>>
    auto sr_data = sr->read_next();
    if (!sr->results_complete()) {
        Rcpp::warning("Read of '%s' incomplete", uri);
    }

    // names gets a vector of strings, at() gets a particular column buffer by name
    auto names = sr_data->get()->names();

    return names;
}

// Initial 'proof of light' function to a column by name
// will likely be replaced / changed in due course
//
//' @noRd
// [[Rcpp::export]]
bool export_column(const std::string& uri, const std::string& colname,
                   SEXP schemaxp, SEXP arrayxp) {
    std::map<std::string, std::string> config;

    // Read all values from the array (return is unique_ptr<SOMAReader>)
    auto sr = tdbs::SOMAReader::open(uri);
    sr->submit();

    // Getting next batch:  std::optinal<std::shared_ptr<ArrayBuffers>>
    auto sr_data = sr->read_next();
    if (!sr->results_complete()) {
        Rcpp::warning("Read of '%s' incomplete", uri);
    }

    // now buf is a shared_ptr to ColumnBuffer
    auto buf = sr_data->get()->at(colname);

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
//' @noRd
// [[Rcpp::export]]
SEXP export_column_direct(const std::string& uri, const std::vector<std::string>& colnames) {

    spdl::info(fmt::format("Reading from {}", uri));

    // Read selected columns from the uri array (return is unique_ptr<SOMAReader>)
    auto sr = tdbs::SOMAReader::open(uri, "", {}, colnames);
    sr->submit();

    // Getting next batch:  std::optinal<std::shared_ptr<ArrayBuffers>>
    auto sr_data = sr->read_next();
    if (!sr->results_complete()) {
        Rcpp::warning("Read of '%s' incomplete", uri);
    }
    spdl::info(fmt::format("Read complete with {} rows and {} cols",
                           sr_data->get()->num_rows(), sr_data->get()->names().size()));

    auto ncol = sr_data->get()->names().size();
    Rcpp::List reslist(ncol);
    reslist.attr("names") = sr_data->get()->names();

    for (size_t i=0; i<ncol; i++) {
        // this allocates, and properly wraps as external pointers controlling lifetime
        SEXP schemaxp = arch_c_allocate_schema();
        SEXP arrayxp = arch_c_allocate_array_data();

        // now buf is a shared_ptr to ColumnBuffer
        auto buf = sr_data->get()->at(colnames[i]);

        spdl::info(fmt::format("Accessing {} at {}", colnames[i], i));

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

//' Read SOMA Data From a Given URI
//'
//' This functions access a given SOMA URI and returns a complete data.frame. It does
//' not iterate; if your data is large than the initial read size consider the \code{sr_*}
//' functions.
//'
//' @param uri Character value with URI path to a SOMA data set
//' @param colnames Optional vector of character value with the name of the columns to retrieve
//' @param qc Optional external Pointer object to TileDB Query Condition, defaults to \sQuote{NULL} i.e.
//' no query condition
//' @param dim_points Optional named list with vector of data points to select on the given
//' dimension(s). Each dimension can be one entry in the list.
//' @param dim_ranges Optional named list with two-column matrix where each row select a range
//' for the given dimension. Each dimension can be one entry in the list.
//' @param loglevel Character value with the desired logging level, defaults to \sQuote{warn}
//' @return An Arrow data structure is returned
//' @examples
//' \dontrun{
//' uri <- "test/soco/pbmc3k_processed/obs"
//' z <- soma_reader(uri)
//' tb <- arrow::as_arrow_table(arch::from_arch_array(z, arrow::RecordBatch))
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List soma_reader(const std::string& uri,
                       Rcpp::Nullable<Rcpp::CharacterVector> colnames = R_NilValue,
                       Rcpp::Nullable<Rcpp::XPtr<tiledb::QueryCondition>> qc = R_NilValue,
                       Rcpp::Nullable<Rcpp::List> dim_points = R_NilValue,
                       Rcpp::Nullable<Rcpp::List> dim_ranges = R_NilValue,
                       const std::string& loglevel = "warn") {

    tdbs::LOG_SET_LEVEL(loglevel);
    spdl::set_level(loglevel);

    spdl::info(fmt::format("[soma_reader] Reading from {}", uri));

    // Read selected columns from the uri (return is unique_ptr<SOMAReader>)
    auto sr = tdbs::SOMAReader::open(uri);

    std::unordered_map<std::string, tiledb_datatype_t> name2type;
    std::shared_ptr<tiledb::ArraySchema> schema = sr->schema();
    tiledb::Domain domain = schema->domain();
    std::vector<tiledb::Dimension> dims = domain.dimensions();
    for (auto& dim: dims) {
        spdl::info(fmt::format("[soma_reader] Dimension {} type {} domain {} extent {}",
                               dim.name(), tiledb::impl::to_str(dim.type()),
                               dim.domain_to_str(), dim.tile_extent_to_str()));
        name2type.emplace(std::make_pair(dim.name(), dim.type()));
    }

    // If we have column names, select them
    if (!colnames.isNull()) {
        std::vector<std::string> cn = Rcpp::as<std::vector<std::string>>(colnames);
        spdl::info(fmt::format("[soma_reader] Selecting {} columns", cn.size()));
        sr->select_columns(cn);
    }

    // If we have a query condition, apply it
    if (!qc.isNull()) {
        spdl::info(fmt::format("[soma_reader] Applying query condition"));
        Rcpp::XPtr<tiledb::QueryCondition> qcxp(qc);
        sr->set_condition(*qcxp);
    }

    // If we have dimension points, apply them
    // The interface is named list, where each (named) list elements is one (named) dimesion
    // The List element is a simple vector of points and each point is applied to the named dimension
    if (!dim_points.isNull()) {
        Rcpp::List lst(dim_points);
        apply_dim_points(sr.get(), name2type, lst);
    }

    // If we have a dimension points, apply them
    if (!dim_ranges.isNull()) {
        Rcpp::List lst(dim_ranges);
        apply_dim_ranges(sr.get(), name2type, lst);
    }

    sr->submit();

    // Getting next batch:  std::optional<std::shared_ptr<ArrayBuffers>>
    auto sr_data = sr->read_next();
    if (!sr->results_complete()) {
        Rcpp::warning("Read of '%s' incomplete", uri);
    }
    spdl::info(fmt::format("[soma_reader] Read complete with {} rows and {} cols",
                               sr_data->get()->num_rows(),
                               sr_data->get()->names().size()));

    const std::vector<std::string> names = sr_data->get()->names();
    auto ncol = names.size();
    Rcpp::List schlst(ncol), arrlst(ncol);

    for (size_t i=0; i<ncol; i++) {
        // this allocates, and properly wraps as external pointers controlling lifetime
        SEXP schemaxp = arch_c_allocate_schema();
        SEXP arrayxp = arch_c_allocate_array_data();

        spdl::info(fmt::format("[soma_reader] Accessing {} at {}", names[i], i));

        // now buf is a shared_ptr to ColumnBuffer
        auto buf = sr_data->get()->at(names[i]);

        // this is pair of array and schema pointer
        auto pp = tdbs::ArrowAdapter::to_arrow(buf);

        memcpy((void*) R_ExternalPtrAddr(schemaxp), pp.second.get(), sizeof(ArrowSchema));
        memcpy((void*) R_ExternalPtrAddr(arrayxp), pp.first.get(), sizeof(ArrowArray));

        spdl::info(fmt::format("[soma_reader] Incoming name {}", std::string(pp.second->name)));

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

//' @noRd
// [[Rcpp::export]]
Rcpp::CharacterVector get_column_types(const std::string& uri,
                                       const std::vector<std::string>& colnames) {

    auto sr = tdbs::SOMAReader::open(uri);
    sr->submit();
    auto sr_data = sr->read_next();
    size_t n = colnames.size();
    Rcpp::CharacterVector vs(n);
    for (size_t i=0; i<n; i++) {
        auto datatype = sr_data->get()->at(colnames[i])->type();
        vs[i] = std::string(tiledb::impl::to_str(datatype));
    }
    vs.attr("names") = colnames;
    return vs;
}

//' @rdname soma_reader
//' @export
// [[Rcpp::export]]
double nnz(const std::string& uri) {
    auto sr = tdbs::SOMAReader::open(uri);
    return static_cast<double>(sr->nnz());
}
