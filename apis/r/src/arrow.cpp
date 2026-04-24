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
void createSchemaForDataFrame(
    const std::string& uri,
    naxpSchema nasp,
    Rcpp::CharacterVector index_column_names,
    Rcpp::List index_column_domains,
    Rcpp::List pclst,
    Rcpp::XPtr<somactx_wrap_t> ctxxp,
    Rcpp::Nullable<Rcpp::DatetimeVector> tsvec) {
    nanoarrow::UniqueSchema sp{nanoarrow_schema_from_xptr(nasp)};
    auto schema = tdbs::common::arrow::make_managed_unique<ArrowSchema>();
    sp.move(schema.get());

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
    pltcfg.override_naming_restriction = Rcpp::as<bool>(pclst["override_naming_restriction"]);

    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // shared pointer to TileDB Context from SOMAContext
    std::shared_ptr<tiledb::Context> ctx = sctx->tiledb_ctx();

    // optional timestamp range
    std::optional<tdbs::TimestampRange> tsrng = makeTimestampRange(tsvec);

    auto encode_domain = [&](std::string_view datatype, Rcpp::RObject domain) -> tdbs::DomainRange {
        auto encode = [&]<typename DomainType, typename ElementType>() -> tdbs::DomainRange {
            if (domain.isNULL()) {
                return std::optional<std::pair<DomainType, DomainType>>();
            }

            auto dom = Rcpp::as<std::vector<ElementType>>(domain);

            if constexpr (std::is_integral_v<DomainType> || std::is_floating_point_v<DomainType>) {
                if (domain.isObject()) {
                    return std::make_optional<std::pair<DomainType, DomainType>>(std::make_pair<DomainType, DomainType>(
                        static_cast<DomainType>(*reinterpret_cast<int64_t*>(&dom[0])),
                        static_cast<DomainType>(*reinterpret_cast<int64_t*>(&dom[1]))));
                } else {
                    return std::make_optional<std::pair<DomainType, DomainType>>(std::make_pair<DomainType, DomainType>(
                        static_cast<DomainType>(dom[0]), static_cast<DomainType>(dom[1])));
                }
            } else {
                return std::make_optional<std::pair<DomainType, DomainType>>(std::make_pair<DomainType, DomainType>(
                    static_cast<DomainType>(dom[0]), static_cast<DomainType>(dom[1])));
            }
        };

        switch (tdbs::common::arrow::to_tiledb_format(datatype)) {
            case TILEDB_CHAR:
            case TILEDB_STRING_ASCII:
            case TILEDB_STRING_UTF8:
            case TILEDB_BLOB:
            case TILEDB_GEOM_WKT:
            case TILEDB_GEOM_WKB:
                return encode.template operator()<std::string, std::string>();
            case TILEDB_INT8:
                return encode.template operator()<int8_t, double_t>();
            case TILEDB_UINT8:
                return encode.template operator()<uint8_t, double_t>();
            case TILEDB_INT16:
                return encode.template operator()<int16_t, double_t>();
            case TILEDB_UINT16:
                return encode.template operator()<uint16_t, double_t>();
            case TILEDB_INT32:
                return encode.template operator()<int32_t, double_t>();
            case TILEDB_UINT32:
                return encode.template operator()<uint32_t, double_t>();
            case TILEDB_DATETIME_SEC:
            case TILEDB_DATETIME_MS:
            case TILEDB_DATETIME_US:
            case TILEDB_DATETIME_NS:
            case TILEDB_INT64:
                return encode.template operator()<int64_t, double_t>();
            case TILEDB_UINT64:
                return encode.template operator()<uint64_t, double_t>();
            case TILEDB_FLOAT32:
                return encode.template operator()<float_t, double_t>();
            case TILEDB_FLOAT64:
                return encode.template operator()<double_t, double_t>();
            default:
                throw std::runtime_error(
                    "[encode_domain] Unsupported type " +
                    tiledb::impl::type_to_str(tdbs::common::arrow::to_tiledb_format(datatype)));
        }
    };

    std::vector<tdbs::DomainRange> domains;
    for (size_t i = 0; i < index_column_names.size(); ++i) {
        std::string column = Rcpp::as<std::string>(index_column_names[i]);
        if (column == tdbs::SOMA_JOINID) {
            domains.emplace_back(
                encode_domain(tdbs::common::arrow::to_arrow_format(TILEDB_INT64), index_column_domains[i]));
            continue;
        }

        bool column_found = false;
        for (int64_t j = 0; j < schema->n_children; ++j) {
            if (schema->children[j]->name == column) {
                domains.emplace_back(encode_domain(schema->children[j]->format, index_column_domains[i]));
                column_found = true;
                break;
            }
        }

        if (!column_found) {
            // If the index column is missing from schema just append an empty pair as domain.
            // The C++ layer is responsible for catching the missing index column from schema and report back the proper error.
            domains.emplace_back(tdbs::DomainRange());
        }
    }

    tdbs::SOMADataFrame::create(
        uri, std::move(schema), Rcpp::as<std::vector<std::string>>(index_column_names), domains, sctx, pltcfg, tsrng);
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
    Rcpp::XPtr<somaobj_wrap_t> soma_array, naxpArray naap, naxpSchema nasp, const std::string arraytype = "") {
    if (!soma_array) {
        Rcpp::stop("Internal error: SOMAObject handle is not initialized.");
    }
    // Move unique arrow array from R to SOMA managed.
    nanoarrow::UniqueArray ap{nanoarrow_array_from_xptr(naap)};
    auto array = tdbs::common::arrow::make_managed_unique<ArrowArray>();
    ap.move(array.get());

    // Move unique arrow schema from R to SOMA managed.
    nanoarrow::UniqueSchema sp{nanoarrow_schema_from_xptr(nasp)};
    auto schema = tdbs::common::arrow::make_managed_unique<ArrowSchema>();
    sp.move(schema.get());

    auto mq = soma_array->ptr<tiledbsoma::SOMAArray>()->create_managed_query("unnamed");
    mq.set_layout(
        arraytype == "SOMADenseNDArray" ? tdbs::common::ResultOrder::colmajor : tdbs::common::ResultOrder::automatic);
    mq.set_array_data(schema.get(), array.get());
    mq.submit_write();
}
