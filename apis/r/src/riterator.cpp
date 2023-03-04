// Work in progress

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
#include "xptr-utils.h"         // xptr taggging utilities

// in rinterface.cpp
SEXP schema_owning_xptr(void);
SEXP array_owning_xptr(void);

namespace tdbs = tiledbsoma;

//' Iterator-Style Access to SOMA Array via SOMAReader
//'
//' The `sr_*` functions provide low-level access to an instance of the SOMAReader
//' class so that iterative access over parts of a (large) array is possible.
//' \describe{
//'   \item{\code{sr_setup}}{instantiates and by default also submits a query}
//'   \item{\code{sr_complete}}{checks if more data is available}
//'   \item{\code{sr_next}}{returns the next chunk}
//' }
//'
//' @param ctx An external pointer to a TileDB Context object
//' @param uri Character value with URI path to a SOMA data set
//' @param colnames Optional vector of character value with the name of the columns to retrieve
//' @param qc Optional external Pointer object to TileDB Query Condition, defaults to \sQuote{NULL} i.e.
//' no query condition
//' @param dim_points Optional named list with vector of data points to select on the given
//' dimension(s). Each dimension can be one entry in the list.
//' @param dim_ranges Optional named list with two-column matrix where each row select a range
//' for the given dimension. Each dimension can be one entry in the list.
//' @param config Optional named chracter vector with \sQuote{key} and \sQuote{value} pairs
//' used as TileDB config parameters. If unset default configuration is used.
//' @param loglevel Character value with the desired logging level, defaults to \sQuote{auto}
//' which lets prior setting prevail, any other value is set as new logging level.
//' @param sr An external pointer to a TileDB SOMAReader object
//'
//' @return \code{sr_setup} returns an external pointer to a SOMAReader. \code{sr_complete}
//' returns a boolean, and \code{sr_next} returns an Arrow array helper object.
//'
//' @examples
//' \dontrun{
//' ctx <- tiledb_ctx()
//' uri <- "test/soco/pbmc3k_processed/obs"
//' sr <- sr_setup(ctx@ptr, uri, "warn")
//' rl <- data.frame()
//' while (nrow(rl) == 0 || !tiledbsoma:::sr_complete(sr)) {
//'     dat <- tiledbsoma:::sr_next(sr)
//'     dat |>
//'         arch::from_arch_array(arrow::RecordBatch) |>
//'         arrow::as_arrow_table() |>
//'         collect() |>
//'         as.data.frame() |>
//'         data.table() -> D
//'     rl <- rbind(rl, D)
//' }
//' summary(rl)
//' }
//' @export
// [[Rcpp::export]]
Rcpp::XPtr<tdbs::SOMAReader> sr_setup(Rcpp::XPtr<tiledb::Context> ctx,
                                      const std::string& uri,
                                      Rcpp::Nullable<Rcpp::CharacterVector> colnames = R_NilValue,
                                      Rcpp::Nullable<Rcpp::XPtr<tiledb::QueryCondition>> qc = R_NilValue,
                                      Rcpp::Nullable<Rcpp::List> dim_points = R_NilValue,
                                      Rcpp::Nullable<Rcpp::List> dim_ranges = R_NilValue,
                                      Rcpp::Nullable<Rcpp::CharacterVector> config = R_NilValue,
                                      const std::string& loglevel = "auto") {
    check_xptr_tag<tiledb::Context>(ctx);
    if (loglevel != "auto") {
        spdl::set_level(loglevel);
        tdbs::LOG_SET_LEVEL(loglevel);
    }

    spdl::info("[sr_setup] Setting up {}", uri);

    std::string_view name = "unnamed";
    std::vector<std::string> column_names = {};
    std::string_view batch_size = "auto";
    std::string_view result_order = "auto";

    std::shared_ptr<tiledb::Context> ctxptr = nullptr;

    std::map<std::string, std::string> platform_config = {};
    if (!config.isNull()) {
        Rcpp::CharacterVector confvec(config);
        Rcpp::CharacterVector namesvec = confvec.attr("names"); // extract names from named R vector
        size_t n = confvec.length();
        for (size_t i = 0; i<n; i++) {
            platform_config.emplace(std::make_pair(std::string(namesvec[i]), std::string(confvec[i])));
            spdl::debug("[sr_setup] config map adding '{}' = '{}'", std::string(namesvec[i]), std::string(confvec[i]));
        }
        tiledb::Config cfg(platform_config);
        spdl::debug("[sr_setup] creating ctx object with supplied config");
        ctxptr = std::make_shared<tiledb::Context>(cfg);
    } else {
        tiledb::Config cfg{ctx.get()->config()}; // get default config in order to make shared_ptr
        spdl::debug("[sr_setup] creating ctx object with default config");
        ctxptr = std::make_shared<tiledb::Context>(cfg);
    }

    if (!colnames.isNull()) {
        column_names = Rcpp::as<std::vector<std::string>>(colnames);
    }

    auto ptr = new tdbs::SOMAReader(uri, name, ctxptr, column_names, batch_size, result_order);

    std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>> name2dim;
    std::shared_ptr<tiledb::ArraySchema> schema = ptr->schema();
    tiledb::Domain domain = schema->domain();
    std::vector<tiledb::Dimension> dims = domain.dimensions();
    for (auto& dim: dims) {
        spdl::info("[soma_reader] Dimension {} type {} domain {} extent {}",
                   dim.name(), tiledb::impl::to_str(dim.type()),
                   dim.domain_to_str(), dim.tile_extent_to_str());
        name2dim.emplace(std::make_pair(dim.name(), std::make_shared<tiledb::Dimension>(dim)));
    }

    // If we have a query condition, apply it
    if (!qc.isNull()) {
        spdl::info("[soma_reader] Applying query condition") ;
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
    Rcpp::XPtr<tdbs::SOMAReader> xptr = make_xptr<tdbs::SOMAReader>(ptr);
    return xptr;
}
//' @rdname sr_setup
//' @export
// [[Rcpp::export]]
bool sr_complete(Rcpp::XPtr<tdbs::SOMAReader> sr) {
   check_xptr_tag<tdbs::SOMAReader>(sr);
   spdl::info("[sr_complete] Complete test is {}", sr->is_complete());
   return sr->is_complete();
}

//' @rdname sr_setup
//' @export
// [[Rcpp::export]]
Rcpp::List sr_next(Rcpp::XPtr<tdbs::SOMAReader> sr) {
   check_xptr_tag<tdbs::SOMAReader>(sr);

   auto sr_data = sr->read_next();
   if (!sr->results_complete()) {
       spdl::trace("[sr_next] Read is incomplete");
   }
   spdl::info("[sr_next] Read {} rows and {} cols",
              sr_data->get()->num_rows(), sr_data->get()->names().size()) ;

   const std::vector<std::string> names = sr_data->get()->names();
   auto ncol = names.size();
   //Rcpp::List schlst(ncol), arrlst(ncol);
   SEXP schemaxp = schema_owning_xptr();
   SEXP arrayxp = array_owning_xptr();
   ArrowSchemaInitFromType((ArrowSchema*)R_ExternalPtrAddr(schemaxp), NANOARROW_TYPE_STRUCT);
   ArrowSchemaAllocateChildren((ArrowSchema*)R_ExternalPtrAddr(schemaxp), ncol);
   ArrowArrayInitFromType((ArrowArray*)R_ExternalPtrAddr(arrayxp), NANOARROW_TYPE_STRUCT);
   ArrowArrayAllocateChildren((ArrowArray*)R_ExternalPtrAddr(arrayxp), ncol);
   
   int data_rows = 0;

   for (size_t i=0; i<ncol; i++) {
       // this allocates, and properly wraps as external pointers controlling lifetime
       SEXP chldschemaxp = schema_owning_xptr();
       SEXP chldarrayxp = array_owning_xptr();

       spdl::trace("[sr_next] Accessing {} at {}", names[i], i);

       // now buf is a shared_ptr to ColumnBuffer
       auto buf = sr_data->get()->at(names[i]);

       // this is pair of array and schema pointer
       auto pp = tdbs::ArrowAdapter::to_arrow(buf);

       memcpy((void*) R_ExternalPtrAddr(chldschemaxp), pp.second.get(), sizeof(ArrowSchema));
       memcpy((void*) R_ExternalPtrAddr(chldarrayxp), pp.first.get(), sizeof(ArrowArray));

       //schlst[i] = schemaxp;
       //arrlst[i] = arrayxp;
       ((ArrowSchema*)R_ExternalPtrAddr(schemaxp))->children[i] = (ArrowSchema*)R_ExternalPtrAddr(chldschemaxp);
       ((ArrowArray*)R_ExternalPtrAddr(arrayxp))->children[i] = (ArrowArray*)R_ExternalPtrAddr(chldarrayxp);

       if (pp.first->length > data_rows) data_rows = pp.first->length;

   }

   ((ArrowArray*)R_ExternalPtrAddr(arrayxp))->length = data_rows;
   spdl::info("[sr_next] Exporting chunk with {} rows", data_rows);
   Rcpp::List as = Rcpp::List::create(Rcpp::Named("array_data") = arrayxp,
                                      Rcpp::Named("schema") = schemaxp);
                                       
   return as;
   
   // struct ArrowArray* array_data_tmp = (struct ArrowArray*) R_ExternalPtrAddr(arrlst[0]);
   // int rows = static_cast<int>(array_data_tmp->length);
   // SEXP sxp = arch_c_schema_xptr_new(Rcpp::wrap("+s"),  // format
   //                                   Rcpp::wrap(""),    // name
   //                                   Rcpp::List(),      // metadata
   //                                   Rcpp::wrap(2),     // flags, 2 == unordered, nullable, no sorted map keys
   //                                   schlst,            // children
   //                                   R_NilValue);       // dictionary
   // SEXP axp = arch_c_array_from_sexp(Rcpp::List::create(Rcpp::Named("")=R_NilValue), // buffers
   //                                   Rcpp::wrap(rows),  // length
   //                                   Rcpp::wrap(-1),    // null count, -1 means not determined
   //                                   Rcpp::wrap(0),     // offset (in bytes)
   //                                   arrlst,            // children
   //                                   R_NilValue);       // dictionary
   // Rcpp::List as = Rcpp::List::create(Rcpp::Named("schema") = sxp,
   //                                    Rcpp::Named("array_data") = axp);
   // as.attr("class") = "arch_array";
   // return as;
}
