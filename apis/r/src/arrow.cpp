#include <Rcpp/Lighter>  // for R interface to C++

#include <nanoarrow/r.h>            // for C/C++ interface to Arrow (via header exported from the R package)
#include <RcppInt64>                // for fromInteger64
#include <nanoarrow/nanoarrow.hpp>  // for C/C++ interface to Arrow (vendored)

#include <tiledbsoma/reindexer/reindexer.h>
#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"  // local declarations
#include "xptr-utils.h"  // xptr taggging utilities

namespace tdbs = tiledbsoma;

void _show_content(const nanoarrow::UniqueArray& ap, const nanoarrow::UniqueSchema& sp) {
    int n = sp.get()->n_children;
    ArrowError ec;
    for (auto i = 0; i < n; i++) {
        Rcpp::Rcout << " " << sp.get()->children[i]->name << " : ";

        nanoarrow::UniqueArrayView col;
        ArrowArrayViewInitFromSchema(col.get(), sp.get()->children[i], &ec);
        ArrowArrayViewSetArray(col.get(), ap.get()->children[i], &ec);

        int m = col.get()->length;
        for (auto j = 0; j < m; j++)
            Rcpp::Rcout << ArrowArrayViewGetDoubleUnsafe(col.get(), j) << " ";
        Rcpp::Rcout << std::endl;
    }
}

//  ctx_wrap_t* ctxwrap_p = new ContextWrapper(ctxptr);
//  Rcpp::XPtr<ctx_wrap_t> ctx_wrap_xptr = make_xptr<ctx_wrap_t>(ctxwrap_p,
//  false);

// [[Rcpp::export]]
void createSchemaFromArrow(
    const std::string& uri,
    naxpSchema nasp,
    naxpArray nadimap,
    naxpSchema nadimsp,
    bool sparse,
    std::string datatype,
    Rcpp::List pclst,
    Rcpp::XPtr<somactx_wrap_t> ctxxp,
    Rcpp::Nullable<Rcpp::DatetimeVector> tsvec = R_NilValue) {
    // struct ArrowArray* ap = (struct ArrowArray*) R_ExternalPtrAddr(naap);
    // struct ArrowSchema* sp = (struct ArrowSchema*) R_ExternalPtrAddr(nasp);
    //
    //  or
    //  auto ap = nanoarrow_array_from_xptr(naap);
    //  auto sp = nanoarrow_schema_from_xptr(nasp);
    //
    //  or:
    // nanoarrow::UniqueArray ap{nanoarrow_array_from_xptr(naap)};
    nanoarrow::UniqueSchema sp{nanoarrow_schema_from_xptr(nasp)};
    //_show_content(ap, sp);
    nanoarrow::UniqueArray apdim{nanoarrow_array_from_xptr(nadimap)};
    nanoarrow::UniqueSchema spdim{nanoarrow_schema_from_xptr(nadimsp)};
    //_show_content(apdim, spdim);

    auto schema = tdbs::common::arrow::make_managed_unique<ArrowSchema>();
    sp.move(schema.get());

    auto dimsch = tdbs::common::arrow::make_managed_unique<ArrowSchema>();
    spdim.move(dimsch.get());
    auto dimarr = tdbs::common::arrow::make_managed_unique<ArrowArray>();
    apdim.move(dimarr.get());

    tdbs::PlatformConfig pltcfg;
    pltcfg.dataframe_dim_zstd_level = Rcpp::as<int>(pclst["dataframe_dim_zstd_level"]);
    pltcfg.sparse_nd_array_dim_zstd_level = Rcpp::as<int>(pclst["sparse_nd_array_dim_zstd_level"]);
    pltcfg.dense_nd_array_dim_zstd_level = Rcpp::as<int>(pclst["dense_nd_array_dim_zstd_level"]);
    pltcfg.write_X_chunked = Rcpp::as<bool>(pclst["write_X_chunked"]);
    pltcfg.goal_chunk_nnz = Rcpp::as<double>(pclst["goal_chunk_nnz"]);
    pltcfg.capacity = Rcpp::as<double>(pclst["capacity"]);
    pltcfg.offsets_filters = Rcpp::as<std::string>(pclst["offsets_filters"]);
    pltcfg.validity_filters = Rcpp::as<std::string>(pclst["validity_filters"]);
    pltcfg.allows_duplicates = Rcpp::as<bool>(pclst["allows_duplicates"]);
    pltcfg.cell_order = Rcpp::as<std::string>(pclst["cell_order"]);
    pltcfg.tile_order = Rcpp::as<std::string>(pclst["tile_order"]);
    pltcfg.attrs = Rcpp::as<std::string>(pclst["attrs"]);
    pltcfg.dims = Rcpp::as<std::string>(pclst["dims"]);

    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // shared pointer to TileDB Context from SOMAContext
    std::shared_ptr<tiledb::Context> ctx = sctx->tiledb_ctx();

    // optional timestamp range
    std::optional<tdbs::TimestampRange> tsrng = makeTimestampRange(tsvec);

    bool exists = false;
    if (datatype == "SOMADataFrame") {
        exists = tdbs::SOMADataFrame::exists(uri, sctx);
    } else if (datatype == "SOMASparseNDArray") {
        exists = tdbs::SOMASparseNDArray::exists(uri, sctx);
    } else if (datatype == "SOMADenseNDArray") {
        exists = tdbs::SOMADenseNDArray::exists(uri, sctx);
    } else {
        Rcpp::stop(tfm::format("Error: Invalid SOMA type_argument '%s'", datatype));
    }

    if (exists) {
        Rcpp::stop(tfm::format("Error: Array '%s' already exists", uri));
    }

    if (datatype == "SOMADataFrame") {
        tdbs::SOMADataFrame::create(
            uri, std::move(schema), std::pair(std::move(dimarr), std::move(dimsch)), sctx, pltcfg, tsrng);
    } else if (datatype == "SOMASparseNDArray") {
        // for arrays n_children will be three as we have two dims and a data
        // col
        std::string datacoltype = sp->children[sp->n_children - 1]->format;
        tdbs::SOMASparseNDArray::create(
            uri, datacoltype, std::pair(std::move(dimarr), std::move(dimsch)), sctx, pltcfg, tsrng);
    } else if (datatype == "SOMADenseNDArray") {
        // for arrays n_children will be three as we have two dims and a data
        // col
        std::string datacoltype = sp->children[sp->n_children - 1]->format;
        tdbs::SOMADenseNDArray::create(
            uri, datacoltype, std::pair(std::move(dimarr), std::move(dimsch)), sctx, pltcfg, tsrng);
    } else {
        Rcpp::stop(tfm::format("Error: Invalid SOMA type_argument '%s'", datatype));
    }
}

// [[Rcpp::export]]
void createSchemaForNDArray(
    const std::string& uri,
    const std::string& format,
    Rcpp::NumericVector shape,
    const std::string& soma_type,
    Rcpp::List pclst,
    Rcpp::XPtr<somactx_wrap_t> ctxxp,
    Rcpp::Nullable<Rcpp::DatetimeVector> tsvec = R_NilValue) {
    tdbs::PlatformConfig pltcfg;
    pltcfg.dataframe_dim_zstd_level = Rcpp::as<int>(pclst["dataframe_dim_zstd_level"]);
    pltcfg.sparse_nd_array_dim_zstd_level = Rcpp::as<int>(pclst["sparse_nd_array_dim_zstd_level"]);
    pltcfg.dense_nd_array_dim_zstd_level = Rcpp::as<int>(pclst["dense_nd_array_dim_zstd_level"]);
    pltcfg.write_X_chunked = Rcpp::as<bool>(pclst["write_X_chunked"]);
    pltcfg.goal_chunk_nnz = Rcpp::as<double>(pclst["goal_chunk_nnz"]);
    pltcfg.capacity = Rcpp::as<double>(pclst["capacity"]);
    pltcfg.offsets_filters = Rcpp::as<std::string>(pclst["offsets_filters"]);
    pltcfg.validity_filters = Rcpp::as<std::string>(pclst["validity_filters"]);
    pltcfg.allows_duplicates = Rcpp::as<bool>(pclst["allows_duplicates"]);
    pltcfg.cell_order = Rcpp::as<std::string>(pclst["cell_order"]);
    pltcfg.tile_order = Rcpp::as<std::string>(pclst["tile_order"]);
    pltcfg.attrs = Rcpp::as<std::string>(pclst["attrs"]);
    pltcfg.dims = Rcpp::as<std::string>(pclst["dims"]);

    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // shared pointer to TileDB Context from SOMAContext
    std::shared_ptr<tiledb::Context> ctx = sctx->tiledb_ctx();

    // optional timestamp range
    std::optional<tdbs::TimestampRange> tsrng = makeTimestampRange(tsvec);

    bool exists = false;
    if (soma_type == "SOMASparseNDArray") {
        exists = tdbs::SOMASparseNDArray::exists(uri, sctx);
    } else if (soma_type == "SOMADenseNDArray") {
        exists = tdbs::SOMADenseNDArray::exists(uri, sctx);
    } else {
        Rcpp::stop(tfm::format("Error: Invalid SOMA type_argument '%s'", soma_type));
    }

    if (exists) {
        Rcpp::stop(tfm::format("Error: Array '%s' already exists", uri));
    }

    std::vector<std::optional<int64_t>> cpp_shape;
    for (size_t i = 0; i < shape.length(); ++i) {
        cpp_shape.push_back(std::make_optional<int64_t>(*reinterpret_cast<int64_t*>(&shape[i])));
    }

    if (soma_type == "SOMASparseNDArray") {
        tdbs::SOMASparseNDArray::create(uri, format, cpp_shape, sctx, pltcfg, tsrng);
    } else if (soma_type == "SOMADenseNDArray") {
        tdbs::SOMADenseNDArray::create(uri, format, cpp_shape, sctx, pltcfg, tsrng);
    } else {
        Rcpp::stop(tfm::format("Error: Invalid SOMA type_argument '%s'", soma_type));
    }
}

// [[Rcpp::export]]
void writeArrayFromArrow(
    Rcpp::XPtr<tiledbsoma::SOMAArray> soma_array, naxpArray naap, naxpSchema nasp, const std::string arraytype = "") {
    if (!soma_array) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    // Move unique arrow array from R to SOMA managed.
    nanoarrow::UniqueArray ap{nanoarrow_array_from_xptr(naap)};
    auto array = tdbs::common::arrow::make_managed_unique<ArrowArray>();
    ap.move(array.get());

    // Move unique arrow schema from R to SOMA managed.
    nanoarrow::UniqueSchema sp{nanoarrow_schema_from_xptr(nasp)};
    auto schema = tdbs::common::arrow::make_managed_unique<ArrowSchema>();
    sp.move(schema.get());

    auto mq = soma_array->create_managed_query("unnamed");
    mq.set_layout(
        arraytype == "SOMADenseNDArray" ? tdbs::common::ResultOrder::colmajor : tdbs::common::ResultOrder::automatic);
    mq.set_array_data(schema.get(), array.get());
    mq.submit_write();
}
