// we currently get deprecation warnings by default which are noisy
#ifndef TILEDB_NO_API_DEPRECATION_WARNINGS
#define TILEDB_NO_API_DEPRECATION_WARNINGS
#endif

#include <Rcpp.h>               // for R interface to C++
#include <nanoarrow.h>          // for C interface to Arrow

#include <tiledb/tiledb>
#if TILEDB_VERSION_MAJOR == 2 && TILEDB_VERSION_MINOR >= 4
#include <tiledb/tiledb_experimental>
#endif

// We get these via nanoarrow and must cannot include carrow.h again
#define ARROW_SCHEMA_AND_ARRAY_DEFINED 1
#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"         // local declarations
#include "xptr-utils.h"         // xptr taggging utilitie
Rcpp::XPtr<ArrowSchema> schema_setup_struct(Rcpp::XPtr<ArrowSchema> schxp, int64_t n_children);
Rcpp::XPtr<ArrowArray> array_setup_struct(Rcpp::XPtr<ArrowArray> arrxp, int64_t n_children);

namespace tdbs = tiledbsoma;

//' Iterator-Style Access to SOMA Array via SOMAArray
//'
//' The `sr_*` functions provide low-level access to an instance of the SOMAArray
//' class so that iterative access over parts of a (large) array is possible.
//' \describe{
//'   \item{\code{sr_setup}}{instantiates and by default also submits a query}
//'   \item{\code{sr_complete}}{checks if more data is available}
//'   \item{\code{sr_next}}{returns the next chunk}
//' }
//'
//' @param uri Character value with URI path to a SOMA data set
//' @param config Named chracter vector with \sQuote{key} and \sQuote{value} pairs
//' used as TileDB config parameters.
//' @param colnames Optional vector of character value with the name of the columns to retrieve
//' @param qc Optional external Pointer object to TileDB Query Condition, defaults to \sQuote{NULL} i.e.
//' no query condition
//' @param dim_points Optional named list with vector of data points to select on the given
//' dimension(s). Each dimension can be one entry in the list.
//' @param dim_ranges Optional named list with two-column matrix where each row select a range
//' for the given dimension. Each dimension can be one entry in the list.
//' @param batch_size Optional argument for size of data batches, defaults to \sQuote{auto}
//' @param result_order Optional argument for query result order, defaults to \sQuote{auto}
//' @param loglevel Character value with the desired logging level, defaults to \sQuote{auto}
//' which lets prior setting prevail, any other value is set as new logging level.
//' @param timestamp_end Optional POSIXct (i.e. Datetime) type for end of interval for which
//' data is considered.
//' @param sr An external pointer to a TileDB SOMAArray object
//'
//' @return \code{sr_setup} returns an external pointer to a SOMAArray. \code{sr_complete}
//' returns a boolean, and \code{sr_next} returns an Arrow array helper object.
//'
//' @examples
//' \dontrun{
//' uri <- extract_dataset("soma-dataframe-pbmc3k-processed-obs")
//' ctx <- tiledb::tiledb_ctx()
//' sr <- sr_setup(uri, config=as.character(tiledb::config(ctx)))
//' rl <- data.frame()
//' while (!sr_complete(sr)) {
//'   dat <- sr_next(sr)
//'   rb <- arrow::RecordBatch$import_from_c(dat$array_data, dat$schema)
//'   rl <- rbind(rl, as.data.frame(rb))
//' }
//' summary(rl)
//' }
//' @export
// [[Rcpp::export]]
Rcpp::XPtr<tdbs::SOMAArray> sr_setup(const std::string& uri,
                                     Rcpp::CharacterVector config,
                                     Rcpp::Nullable<Rcpp::CharacterVector> colnames = R_NilValue,
                                     Rcpp::Nullable<Rcpp::XPtr<tiledb::QueryCondition>> qc = R_NilValue,
                                     Rcpp::Nullable<Rcpp::List> dim_points = R_NilValue,
                                     Rcpp::Nullable<Rcpp::List> dim_ranges = R_NilValue,
                                     std::string batch_size = "auto",
                                     std::string result_order = "auto",
                                     Rcpp::Nullable<Rcpp::Datetime> timestamp_end = R_NilValue,
                                     const std::string& loglevel = "auto") {

    if (loglevel != "auto") {
        spdl::set_level(loglevel);
        tdbs::LOG_SET_LEVEL(loglevel);
    }

    spdl::debug("[sr_setup] Setting up {}", uri);

    std::string_view name = "unnamed";
    std::vector<std::string> column_names = {};

    std::shared_ptr<tiledb::Context> ctxptr = nullptr;

    std::map<std::string, std::string> platform_config = config_vector_to_map(Rcpp::wrap(config));
    tiledb::Config cfg(platform_config);
    spdl::debug("[sr_setup] creating ctx object with supplied config");
    ctxptr = std::make_shared<tiledb::Context>(cfg);

    if (!colnames.isNull()) {
        column_names = Rcpp::as<std::vector<std::string>>(colnames);
    }

    std::uint64_t ts_start = 0; 		// beginning of time aka the epoch (force double signature)
    std::uint64_t ts_end = std::numeric_limits<uint64_t>::max(); 	// max if unset
    if (!timestamp_end.isNull()) {
        ts_end = Rcpp::as<Rcpp::Datetime>(timestamp_end).getFractionalTimestamp() * 1e3; // in msec
        spdl::info(tfm::format("[sr_setup] ts_end set to %ld", ts_end));
    }

    auto tdb_result_order = get_tdb_result_order(result_order);

    auto ptr = new tdbs::SOMAArray(OpenMode::read, uri, name, ctxptr, column_names, batch_size,
                                   tdb_result_order, std::make_pair(ts_start, ts_end));

    std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>> name2dim;
    std::shared_ptr<tiledb::ArraySchema> schema = ptr->schema();
    tiledb::Domain domain = schema->domain();
    std::vector<tiledb::Dimension> dims = domain.dimensions();
    for (auto& dim: dims) {
        spdl::debug("[soma_array_reader] Dimension {} type {} domain {} extent {}",
                    dim.name(), tiledb::impl::to_str(dim.type()),
                    dim.domain_to_str(), dim.tile_extent_to_str());
        name2dim.emplace(std::make_pair(dim.name(), std::make_shared<tiledb::Dimension>(dim)));
    }

    // If we have a query condition, apply it
    if (!qc.isNull()) {
        spdl::debug("[soma_array_reader] Applying query condition") ;
        Rcpp::XPtr<tiledb::QueryCondition> qcxp(qc);
        ptr->set_condition(*qcxp);
    }

    // If we have dimension points, apply them
    // The interface is named list, where each (named) list elements is one (named) dimesion
    // The List element is a simple vector of points and each point is applied to the named dimension
    if (!dim_points.isNull()) {
        Rcpp::List lst(dim_points);
        apply_dim_points(ptr, name2dim, lst);
    }

    // If we have a dimension points, apply them
    if (!dim_ranges.isNull()) {
        Rcpp::List lst(dim_ranges);
        apply_dim_ranges(ptr, name2dim, lst);
    }

    ptr->submit();
    Rcpp::XPtr<tdbs::SOMAArray> xptr = make_xptr<tdbs::SOMAArray>(ptr);
    return xptr;
}

//' @rdname sr_setup
//' @export
// [[Rcpp::export]]
bool sr_complete(Rcpp::XPtr<tdbs::SOMAArray> sr) {
   check_xptr_tag<tdbs::SOMAArray>(sr);
   bool complt = sr->is_complete(true);
   bool initial = sr->is_initial_read();
   bool res = complt && !initial; // completed transfer if query status complete and query ran once
   spdl::debug("[sr_complete] Complete query test {} (compl {} initial {})", res, complt, initial);
   return res;
}

Rcpp::List create_empty_arrow_table() {
    Rcpp::XPtr<ArrowSchema> schemaxp = schema_owning_xptr();
    Rcpp::XPtr<ArrowArray> arrayxp = array_owning_xptr();
    schemaxp = schema_setup_struct(schemaxp, 0);
    arrayxp = array_setup_struct(arrayxp, 0);
    arrayxp->length = 0;
    Rcpp::List as = Rcpp::List::create(Rcpp::Named("array_data") = arrayxp,
                                       Rcpp::Named("schema") = schemaxp);
    return as;
}


//' @rdname sr_setup
//' @export
// [[Rcpp::export]]
Rcpp::List sr_next(Rcpp::XPtr<tdbs::SOMAArray> sr) {
   check_xptr_tag<tdbs::SOMAArray>(sr);

   if (sr_complete(sr)) {
       spdl::trace("[sr_next] complete {} num_cells {}",
                   sr->is_complete(true), sr->total_num_cells());
       return create_empty_arrow_table();
   }

   if (!sr->is_initial_read() && sr->total_num_cells() == 0) {
       spdl::trace("[sr_next] is_initial_read {} num_cells {}",
                   sr->is_initial_read(), sr->total_num_cells());
       return create_empty_arrow_table();
   }

   auto sr_data = sr->read_next();
   spdl::debug("[sr_next] Read {} rows and {} cols",
               sr_data->get()->num_rows(), sr_data->get()->names().size());

   const std::vector<std::string> names = sr_data->get()->names();
   auto ncol = names.size();
   Rcpp::XPtr<ArrowSchema> schemaxp = schema_owning_xptr();
   Rcpp::XPtr<ArrowArray> arrayxp = array_owning_xptr();
   schemaxp = schema_setup_struct(schemaxp, ncol);
   arrayxp = array_setup_struct(arrayxp, ncol);
   arrayxp->length = 0;

   for (size_t i=0; i<ncol; i++) {
       // this allocates, and properly wraps as external pointers controlling lifetime
       Rcpp::XPtr<ArrowSchema> chldschemaxp = schema_owning_xptr();
       Rcpp::XPtr<ArrowArray> chldarrayxp = array_owning_xptr();

       spdl::trace("[sr_next] Accessing {} at {}", names[i], i);

       // now buf is a shared_ptr to ColumnBuffer
       auto buf = sr_data->get()->at(names[i]);

       // this is pair of array and schema pointer
       auto pp = tdbs::ArrowAdapter::to_arrow(buf);

       memcpy((void*) chldschemaxp, pp.second.get(), sizeof(ArrowSchema));
       memcpy((void*) chldarrayxp, pp.first.get(), sizeof(ArrowArray));

       schemaxp->children[i] = chldschemaxp;
       arrayxp->children[i] = chldarrayxp;

       if (pp.first->length > arrayxp->length) {
           spdl::debug("[soma_array_reader] Setting array length to {}", pp.first->length);
           arrayxp->length = pp.first->length;
       }

   }

   spdl::debug("[sr_next] Exporting chunk with {} rows", arrayxp->length);
   Rcpp::List as = Rcpp::List::create(Rcpp::Named("array_data") = arrayxp,
                                      Rcpp::Named("schema") = schemaxp);
   return as;
}
