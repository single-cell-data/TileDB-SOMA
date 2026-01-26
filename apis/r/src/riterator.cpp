// we currently get deprecation warnings by default which are noisy
#ifndef TILEDB_NO_API_DEPRECATION_WARNINGS
#define TILEDB_NO_API_DEPRECATION_WARNINGS
#endif

#include <Rcpp.h>  // for R interface to C++
#include <nanoarrow/nanoarrow.h>
#include <nanoarrow/r.h>  // for C interface to Arrow (via R package nanoarrow)
#include <RcppInt64>      // for fromInteger64
#include <sstream>

#include <tiledb/tiledb>
#if TILEDB_VERSION_MAJOR == 2 && TILEDB_VERSION_MINOR >= 4
#include <tiledb/tiledb_experimental>
#endif
#include <tiledbsoma/reindexer/reindexer.h>
#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"  // local declarations
#include "xptr-utils.h"  // xptr taggging utilitie

namespace tdbs = tiledbsoma;

// clang-format off
// Iterator-Style Access to SOMA Array via ManagedQuery
//
// The `mq_*` functions provide low-level access to an instance of a
// ManagedQuery class so that iterative access over parts of a (large) array
// is possible.
// \describe{
//   \item{\code{mq_setup}}{instantiates and by default also submits a query}
//   \item{\code{mq_complete}}{checks if more data is available}
//   \item{\code{mq_next}}{returns the next chunk}
// }
//
// @param uri Character value with URI path to a SOMA data set
// @param config Named chracter vector with \dQuote{key} and \dQuote{value}
// pairs used as TileDB config parameters.
// @param colnames Optional vector of character value with the name of the
// columns to retrieve
// @param qc Optional external Pointer object to TileDB Query Condition,
// defaults to \code{NULL} i.e. no query condition
// @param dim_points Optional named list with vector of data points to select
// on the given dimension(s). Each dimension can be one entry in the list.
// @param dim_ranges Optional named list with two-column matrix where each row
// select a range for the given dimension. Each dimension can be one entry in
// the list.
// @param batch_size Optional argument for size of data batches, defaults to
// \dQuote{\code{auto}}
// @param result_order Optional argument for query result order, defaults to
// \dQuote{\code{auto}}
// @param loglevel Character value with the desired logging level, defaults to
// \dQuote{\code{auto}} which lets prior setting prevail, any other value is
// set as new logging level.
// @param timestamprange Optional POSIXct (i.e. Datetime) vector with start
// and end of interval for which data is considered.
// @param mq An external pointer to a ManagedQuery object.
//
// @return \code{mq_setup}: returns an external pointer to a ManagedQuery
// object.
//
// @examples
// \dontrun{
// uri <- extract_dataset("soma-dataframe-pbmc3k-processed-obs")
// ctxcp <- soma_context()
// mq <- mq_setup(uri, ctxxp)
// rl <- data.frame()
// while (!mq_complete(mq)) {
//   dat <- mq_next(mq)
//   rb <- arrow::RecordBatch$import_from_c(dat$array_data, dat$schema)
//   rl <- rbind(rl, as.data.frame(rb))
// }
// summary(rl)
// }
//
// @noRd
// clang-format on
// [[Rcpp::export]]
Rcpp::XPtr<tdbs::ManagedQuery> mq_setup(
    const std::string& uri,
    Rcpp::XPtr<somactx_wrap_t> ctxxp,
    Rcpp::Nullable<Rcpp::CharacterVector> colnames = R_NilValue,
    Rcpp::Nullable<Rcpp::XPtr<tiledb::QueryCondition>> qc = R_NilValue,
    Rcpp::Nullable<Rcpp::List> dim_points = R_NilValue,
    Rcpp::Nullable<Rcpp::List> dim_ranges = R_NilValue,
    std::string batch_size = "auto",
    std::string result_order = "auto",
    Rcpp::Nullable<Rcpp::DatetimeVector> timestamprange = R_NilValue,
    const std::string& loglevel = "auto") {
    if (loglevel != "auto") {
        tdbs::LOG_SET_LEVEL(loglevel);
    }

    std::stringstream ss;
    ss << "[mq_setup] Setting up " << uri;
    tdbs::LOG_DEBUG(ss.str());

    std::string_view name = "unnamed";
    std::vector<std::string> column_names = {};

    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> somactx = ctxxp->ctxptr;

    if (!colnames.isNull()) {
        column_names = Rcpp::as<std::vector<std::string>>(colnames);
    }

    // optional timestamp range
    std::optional<tdbs::TimestampRange> tsrng = makeTimestampRange(timestamprange);

    auto tdb_result_order = get_tdb_result_order(result_order);

    auto arr = tdbs::SOMAArray(OpenMode::soma_read, uri, somactx, tsrng);

    auto mq = new tdbs::ManagedQuery(arr.tiledb_array(), somactx->tiledb_ctx(), name);
    mq->set_layout(tdb_result_order);
    if (!column_names.empty()) {
        mq->select_columns(column_names);
    }

    std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>> name2dim;
    std::shared_ptr<tiledb::ArraySchema> schema = arr.tiledb_schema();
    tiledb::Domain domain = schema->domain();
    std::vector<tiledb::Dimension> dims = domain.dimensions();
    for (auto& dim : dims) {
        std::stringstream ss;
        ss << "[mq_setup] Dimension '" << dim.name() << "' type " << tiledb::impl::to_str(dim.type()) << " domain "
           << dim.domain_to_str() << " extent " << dim.tile_extent_to_str();
        tdbs::LOG_DEBUG(ss.str());
        name2dim.emplace(std::make_pair(dim.name(), std::make_shared<tiledb::Dimension>(dim)));
    }

    // If we have a query condition, apply it
    if (!qc.isNull()) {
        tdbs::LOG_DEBUG("[mq_setup] Applying query condition");
        Rcpp::XPtr<tiledb::QueryCondition> qcxp(qc);
        mq->set_condition(*qcxp);
    }

    // If we have dimension points, apply them
    // The interface is named list, where each (named) list elements is one
    // (named) dimesion The List element is a simple vector of points and each
    // point is applied to the named dimension
    if (!dim_points.isNull()) {
        Rcpp::List lst(dim_points);
        apply_dim_points(mq, name2dim, lst);
    }

    // If we have a dimension points, apply them
    if (!dim_ranges.isNull()) {
        Rcpp::List lst(dim_ranges);
        apply_dim_ranges(mq, name2dim, lst);
    }

    Rcpp::XPtr<tdbs::ManagedQuery> xptr = make_xptr<tdbs::ManagedQuery>(mq);
    return xptr;
}

// @rdname mq_setup
//
// @return \code{mq_complete}: returns a boolean
//
// @noRd
//
// [[Rcpp::export]]
bool mq_complete(Rcpp::XPtr<tdbs::ManagedQuery> mq) {
    check_xptr_tag<tdbs::ManagedQuery>(mq);
    return mq->results_complete();
}

// [[Rcpp::export]]
SEXP create_empty_arrow_table() {
    int ncol = 0;

    // Schema first
    auto schemaxp = nanoarrow_schema_owning_xptr();
    auto sch = nanoarrow_output_schema_from_xptr(schemaxp);
    exitIfError(ArrowSchemaInitFromType(sch, NANOARROW_TYPE_STRUCT), "Bad schema init");
    exitIfError(ArrowSchemaSetName(sch, ""), "Bad schema name");
    exitIfError(ArrowSchemaAllocateChildren(sch, ncol), "Bad schema children alloc");

    // Array second
    auto arrayxp = nanoarrow_array_owning_xptr();
    auto arr = nanoarrow_output_array_from_xptr(arrayxp);
    exitIfError(ArrowArrayInitFromType(arr, NANOARROW_TYPE_STRUCT), "Bad array init");
    exitIfError(ArrowArrayAllocateChildren(arr, ncol), "Bad array children alloc");
    arr->length = 0;

    // Nanoarrow special: stick schema into xptr tag to return single SEXP
    array_xptr_set_schema(arrayxp, schemaxp);  // embed schema in array

    return arrayxp;
}

// @rdname mq_setup
//
// @return \code{mq_next}: returns an Arrow array helper object.
//
// @noRd
//
// [[Rcpp::export]]
SEXP mq_next(Rcpp::XPtr<tdbs::ManagedQuery> mq) {
    check_xptr_tag<tdbs::ManagedQuery>(mq);

    if (mq_complete(mq)) {
        std::stringstream ss;
        ss << "[mq_next] complete " << mq->is_complete(true) << " num_cells " << mq->total_num_cells();
        tdbs::LOG_TRACE(ss.str());
        return create_empty_arrow_table();
    }

    auto mq_data = mq->read_next();
    {
        std::stringstream ss;
        ss << "[mq_next] Read " << mq_data->get()->num_rows() << " rows and " << mq_data->get()->names().size()
           << " cols";
        tdbs::LOG_DEBUG(ss.str());
    }

    if (!mq_data) {
        tdbs::LOG_TRACE("[mq_next] complete - mq_data read no data");
        return create_empty_arrow_table();
    }

    const std::vector<std::string> names = mq_data->get()->names();
    auto ncol = names.size();
    // Schema first
    auto schemaxp = nanoarrow_schema_owning_xptr();
    auto sch = nanoarrow_output_schema_from_xptr(schemaxp);
    exitIfError(ArrowSchemaInitFromType(sch, NANOARROW_TYPE_STRUCT), "Bad schema init");
    exitIfError(ArrowSchemaSetName(sch, ""), "Bad schema name");
    exitIfError(ArrowSchemaAllocateChildren(sch, ncol), "Bad schema children alloc");

    // Array second
    auto arrayxp = nanoarrow_array_owning_xptr();
    auto arr = nanoarrow_output_array_from_xptr(arrayxp);
    exitIfError(ArrowArrayInitFromType(arr, NANOARROW_TYPE_STRUCT), "Bad array init");
    exitIfError(ArrowArrayAllocateChildren(arr, ncol), "Bad array children alloc");

    arr->length = 0;  // initial value

    auto pp_vector = tdbs::ArrowAdapter::buffer_to_arrow(mq_data.value(), true);

    for (size_t i = 0; i < ncol; i++) {
        {
            std::stringstream ss;
            ss << "[mq_next] Accessing " << names[i] << " at " << i;
            tdbs::LOG_TRACE(ss.str());
        }

        // this is pair of array and schema pointer
        auto& pp = pp_vector[i];

        ArrowArrayMove(pp.first.get(), arr->children[i]);
        ArrowSchemaMove(pp.second.get(), sch->children[i]);

        if (pp.first->length > arr->length) {
            std::stringstream ss;
            ss << "[soma_array_reader] Setting array length to " << pp.first->length;
            tdbs::LOG_DEBUG(ss.str());
            arr->length = pp.first->length;
        }
    }

    std::stringstream ss;
    ss << "[mq_next] Exporting chunk with " << arr->length << " rows";
    tdbs::LOG_DEBUG(ss.str());
    // Nanoarrow special: stick schema into xptr tag to return single SEXP
    array_xptr_set_schema(arrayxp, schemaxp);  // embed schema in array
    return arrayxp;
}

// [[Rcpp::export]]
void mq_reset(Rcpp::XPtr<tdbs::ManagedQuery> mq) {
    check_xptr_tag<tdbs::ManagedQuery>(mq);
    mq->reset();
    tdbs::LOG_DEBUG("[mq_reset] Reset SOMAArray object");
}

// [[Rcpp::export]]
void mq_set_dim_points(Rcpp::XPtr<tdbs::ManagedQuery> mq, std::string dim, Rcpp::NumericVector points) {
    check_xptr_tag<tdbs::ManagedQuery>(mq);
    // check args ?

    std::vector<int64_t> vec = Rcpp::fromInteger64(points);
    mq->select_points<int64_t>(dim, vec);
    std::stringstream ss;
    ss << "[mq_set_dim_points] Set on dim '" << dim << "' for " << points.length() << " points, first two are "
       << vec[0] << " and " << vec[1];
    tdbs::LOG_DEBUG(ss.str());
}
