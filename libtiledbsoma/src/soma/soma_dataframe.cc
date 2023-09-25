/**
 * @file   soma_dataframe.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2023 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMADataFrame class.
 */

#include <filesystem>

#include <tiledb/tiledb>
#include "array_buffers.h"
#include "soma_dataframe.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

// TileDB type information
struct TypeInfo {
    tiledb_datatype_t type;
    uint64_t elem_size;
    uint32_t cell_val_num;

    // is this represented as "Arrow large"
    bool arrow_large;
};

TypeInfo arrow_type_to_tiledb(ArrowSchema* arw_schema) {
    auto fmt = std::string(arw_schema->format);
    bool large = false;
    if (fmt == "+l") {
        large = false;
        assert(arw_schema->n_children == 1);
        arw_schema = arw_schema->children[0];
    } else if (fmt == "+L") {
        large = true;
        assert(arw_schema->n_children == 1);
        arw_schema = arw_schema->children[0];
    }

    if (fmt == "i")
        return {TILEDB_INT32, 4, 1, large};
    else if (fmt == "l")
        return {TILEDB_INT64, 8, 1, large};
    else if (fmt == "f")
        return {TILEDB_FLOAT32, 4, 1, large};
    else if (fmt == "g")
        return {TILEDB_FLOAT64, 8, 1, large};
    else if (fmt == "b")
        return {TILEDB_BOOL, 1, 1, large};
    else if (fmt == "B")
        return {TILEDB_BLOB, 1, 1, large};
    else if (fmt == "c")
        return {TILEDB_INT8, 1, 1, large};
    else if (fmt == "C")
        return {TILEDB_UINT8, 1, 1, large};
    else if (fmt == "s")
        return {TILEDB_INT16, 2, 1, large};
    else if (fmt == "S")
        return {TILEDB_UINT16, 2, 1, large};
    else if (fmt == "I")
        return {TILEDB_UINT32, 4, 1, large};
    else if (fmt == "L")
        return {TILEDB_UINT64, 8, 1, large};
    // this is kind of a hack
    // technically 'tsn:' is timezone-specific, which we don't support
    // however, the blank (no suffix) base is interconvertible w/
    // np.datetime64
    else if (fmt == "tsn:")
        return {TILEDB_DATETIME_NS, 8, 1, large};
    else if (fmt == "z" || fmt == "Z")
        return {TILEDB_CHAR, 1, TILEDB_VAR_NUM, fmt == "Z"};
    else if (fmt == "u" || fmt == "U")
        return {TILEDB_STRING_UTF8, 1, TILEDB_VAR_NUM, fmt == "U"};
    else
        throw tiledb::TileDBError(
            "[TileDB-Arrow]: Unknown or unsupported Arrow format string '" +
            fmt + "'");
};

template <typename T>
std::vector<std::byte> __aux_schema_get_full_domain(T min, T max) {
    std::byte* min_ptr = reinterpret_cast<std::byte*>(&min);
    std::byte* max_ptr = reinterpret_cast<std::byte*>(&max);
    std::vector<std::byte> slot_domain(min_ptr, min_ptr + sizeof(T));
    slot_domain.insert(slot_domain.end(), max_ptr, max_ptr + sizeof(T));
    return slot_domain;
}

template <typename T>
std::vector<std::byte> __aux_schema_get_extent(T tile) {
    std::byte* tile_ptr = reinterpret_cast<std::byte*>(&tile);
    return std::vector<std::byte>(tile_ptr, tile_ptr + sizeof(T));
}

std::vector<std::byte> _schema_get_full_domain(ArrowSchema* arw_schema) {
    auto datatype = arrow_type_to_tiledb(arw_schema).type;
    switch (datatype) {
            // case TILEDB_STRING_ASCII:
            //     break;

        case TILEDB_FLOAT32: {
            float min = std::numeric_limits<float>::min();
            float max = std::numeric_limits<float>::max();
            return __aux_schema_get_full_domain(min, max);
        }
        case TILEDB_FLOAT64: {
            double min = std::numeric_limits<double>::min();
            double max = std::numeric_limits<double>::max();
            return __aux_schema_get_full_domain(min, max);
        }
        case TILEDB_UINT8: {
            uint8_t min = std::numeric_limits<uint8_t>::min();
            uint8_t max = std::numeric_limits<uint8_t>::max();
            return __aux_schema_get_full_domain(min, max);
        }
        case TILEDB_INT8: {
            int8_t min = std::numeric_limits<int8_t>::min();
            int8_t max = std::numeric_limits<int8_t>::max();
            return __aux_schema_get_full_domain(min, max);
        }
        case TILEDB_UINT16: {
            uint16_t min = std::numeric_limits<uint16_t>::min();
            uint16_t max = std::numeric_limits<uint16_t>::max();
            return __aux_schema_get_full_domain(min, max);
        }
        case TILEDB_INT16: {
            int16_t min = std::numeric_limits<int16_t>::min();
            int16_t max = std::numeric_limits<int16_t>::max();
            return __aux_schema_get_full_domain(min, max);
        }
        case TILEDB_UINT32: {
            uint32_t min = std::numeric_limits<uint32_t>::min();
            uint32_t max = std::numeric_limits<uint32_t>::max();
            return __aux_schema_get_full_domain(min, max);
        }
        case TILEDB_INT32: {
            int32_t min = std::numeric_limits<int32_t>::min();
            int32_t max = std::numeric_limits<int32_t>::max();
            return __aux_schema_get_full_domain(min, max);
        }
        case TILEDB_UINT64: {
            uint64_t min = std::numeric_limits<uint64_t>::min();
            uint64_t max = std::numeric_limits<uint64_t>::max();
            return __aux_schema_get_full_domain(min, max);
        }
        case TILEDB_INT64:
        case TILEDB_DATETIME_YEAR:
        case TILEDB_DATETIME_MONTH:
        case TILEDB_DATETIME_WEEK:
        case TILEDB_DATETIME_DAY:
        case TILEDB_DATETIME_HR:
        case TILEDB_DATETIME_MIN:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_DATETIME_PS:
        case TILEDB_DATETIME_FS:
        case TILEDB_DATETIME_AS: {
            int64_t min, max;
            if (strcmp(arw_schema->name, "soma_joinid") == 0) {
                min = std::numeric_limits<uint32_t>::min();
                max = std::numeric_limits<uint32_t>::max();
            } else {
                min = std::numeric_limits<int64_t>::min();
                max = std::numeric_limits<int64_t>::max();
            }
            return __aux_schema_get_full_domain(min, max);
        }
        default:
            throw tiledb::TileDBError(
                "[TileDB-Arrow]: Unsupported TileDB type for dimension)");
    }
};

std::vector<std::byte> _schema_get_extent(ArrowSchema* arw_schema) {
    auto datatype = arrow_type_to_tiledb(arw_schema).type;
    switch (datatype) {
            // case TILEDB_STRING_ASCII:
            //     break;

        case TILEDB_FLOAT32: {
            float tile = 2048;
            return __aux_schema_get_extent(tile);
        }
        case TILEDB_FLOAT64: {
            double tile = 2048;
            return __aux_schema_get_extent(tile);
        }
        case TILEDB_UINT8: {
            uint8_t tile = 64;
            return __aux_schema_get_extent(tile);
        }
        case TILEDB_INT8: {
            int8_t tile = 64;
            return __aux_schema_get_extent(tile);
        }
        case TILEDB_UINT16: {
            uint16_t tile = 2048;
            return __aux_schema_get_extent(tile);
        }
        case TILEDB_INT16: {
            int16_t tile = 2048;
            return __aux_schema_get_extent(tile);
        }
        case TILEDB_UINT32: {
            uint32_t tile = 2048;
            return __aux_schema_get_extent(tile);
        }
        case TILEDB_INT32: {
            int32_t tile = 2048;
            return __aux_schema_get_extent(tile);
        }
        case TILEDB_UINT64: {
            uint64_t tile = 2048;
            return __aux_schema_get_extent(tile);
        }
        case TILEDB_INT64:
        case TILEDB_DATETIME_YEAR:
        case TILEDB_DATETIME_MONTH:
        case TILEDB_DATETIME_WEEK:
        case TILEDB_DATETIME_DAY:
        case TILEDB_DATETIME_HR:
        case TILEDB_DATETIME_MIN:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_DATETIME_PS:
        case TILEDB_DATETIME_FS:
        case TILEDB_DATETIME_AS: {
            int64_t tile = 2048;
            return __aux_schema_get_extent(tile);
        }
        default:
            throw tiledb::TileDBError(
                "[TileDB-Arrow]: Unsupported TileDB type for dimension)");
    }
};

std::unique_ptr<SOMADataFrame> SOMADataFrame::create(
    std::string_view uri,
    ArrowSchema* schema,
    std::vector<std::string> index_column_names,
    std::map<std::string, std::string> platform_config,
    std::vector<ArrowArray*> domain) {
    auto ctx = std::make_shared<Context>(Config(platform_config));

    if (domain.size() != index_column_names.size()) {
        throw TileDBSOMAError(
            "if domain is specified, it must have the same length as "
            "index_column_names");
    }

    ArraySchema tdb_schema(*ctx, TILEDB_SPARSE);
    Domain tdb_dom(*ctx);

    for (auto index_column : index_column_names) {
        for (int64_t i = 0; i < schema->n_children; ++i) {
            auto child = schema->children[i];
            if (child->name == index_column) {
                auto typeinfo = arrow_type_to_tiledb(child);
                std::vector<std::byte> slot_domain;
                std::vector<std::byte> tile_extent;

                if (domain.size() != 0) {
                    // EXTRACT THIS FROM THE DOMAIN ARROW TABLE
                    auto dom = domain[i];
                    auto d = domain->children[0].buffers[0];
                    // slot_domain = domain->children[i].buffers[0];
                    // slot_domain = nullptr;
                    // tile_extent = nullptr;
                    ;
                } else {
                    slot_domain = _schema_get_full_domain(child);
                    tile_extent = _schema_get_extent(child);
                }

                tdb_dom.add_dimension(Dimension::create(
                    *ctx,
                    child->name,
                    typeinfo.type,
                    slot_domain.data(),
                    tile_extent.data()));
            }
        }
    }
    tdb_schema.set_domain(tdb_dom);

    for (int64_t i = 0; i < schema->n_children; ++i) {
        auto child = schema->children[i];
        if (std::find(
                index_column_names.begin(),
                index_column_names.end(),
                child->name) == index_column_names.end()) {
            auto typeinfo = arrow_type_to_tiledb(child);
            auto attr = Attribute(*ctx, child->name, typeinfo.type);
            tdb_schema.add_attribute(attr);
        }
    }

    // remember enums

    SOMAArray::create(ctx, uri, tdb_schema, "SOMADataFrame");
    return SOMADataFrame::open(uri, OpenMode::read, ctx);
}  // namespace tiledbsoma

std::unique_ptr<SOMADataFrame> SOMADataFrame::create(
    std::string_view uri,
    ArraySchema schema,
    std::map<std::string, std::string> platform_config) {
    return SOMADataFrame::create(
        uri, schema, std::make_shared<Context>(Config(platform_config)));
}

std::unique_ptr<SOMADataFrame> SOMADataFrame::create(
    std::string_view uri, ArraySchema schema, std::shared_ptr<Context> ctx) {
    SOMAArray::create(ctx, uri, schema, "SOMADataFrame");
    return SOMADataFrame::open(uri, OpenMode::read, ctx);
}

std::unique_ptr<SOMADataFrame> SOMADataFrame::open(
    std::string_view uri,
    OpenMode mode,
    std::map<std::string, std::string> platform_config,
    std::vector<std::string> column_names,
    ResultOrder result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    return SOMADataFrame::open(
        uri,
        mode,
        std::make_shared<Context>(Config(platform_config)),
        column_names,
        result_order,
        timestamp);
}

std::unique_ptr<SOMADataFrame> SOMADataFrame::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<Context> ctx,
    std::vector<std::string> column_names,
    ResultOrder result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    return std::make_unique<SOMADataFrame>(
        mode, uri, ctx, column_names, result_order, timestamp);
}

//===================================================================
//= public non-static
//===================================================================

SOMADataFrame::SOMADataFrame(
    OpenMode mode,
    std::string_view uri,
    std::shared_ptr<Context> ctx,
    std::vector<std::string> column_names,
    ResultOrder result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    std::string array_name = std::filesystem::path(uri).filename();
    array_ = std::make_shared<SOMAArray>(
        mode,
        uri,
        array_name,  // label used when debugging
        ctx,
        column_names,
        "auto",  // batch_size,
        result_order,
        timestamp);
    array_->reset();
    array_->submit();
}

void SOMADataFrame::open(
    OpenMode mode, std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    array_->open(mode, timestamp);
    array_->reset();
    array_->submit();
}

void SOMADataFrame::close() {
    array_->close();
}

bool SOMADataFrame::is_open() const {
    return array_->is_open();
}

const std::string SOMADataFrame::uri() const {
    return array_->uri();
}

std::shared_ptr<Context> SOMADataFrame::ctx() {
    return array_->ctx();
}

std::shared_ptr<ArraySchema> SOMADataFrame::schema() const {
    return array_->schema();
}

const std::vector<std::string> SOMADataFrame::index_column_names() const {
    return array_->dimension_names();
}

int64_t SOMADataFrame::count() const {
    return array_->ndim();
}

std::optional<std::shared_ptr<ArrayBuffers>> SOMADataFrame::read_next() {
    return array_->read_next();
}

void SOMADataFrame::write(std::shared_ptr<ArrayBuffers> buffers) {
    array_->reset();
    array_->submit();
    array_->write(buffers);
}

}  // namespace tiledbsoma
