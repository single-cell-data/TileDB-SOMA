/**
 * @file   util.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines utilities.
 */

#include "utils/util.h"
#include <cstring>
#include "logger.h"
#include "utils/logger.h"  // for fmt::format

namespace tiledbsoma::util {

template <typename T>
VarlenBufferPair to_varlen_buffers(std::vector<T> data, bool arrow) {
    size_t nbytes = 0;
    for (auto& elem : data) {
        nbytes += elem.size();
    }

    std::string result;
    std::vector<uint64_t> offsets(data.size() + 1);
    size_t offset = 0;
    size_t idx = 0;

    for (auto& elem : data) {
        result += elem;
        offsets[idx++] = offset;
        offset += elem.size();
    }
    offsets[idx] = offset;

    // Remove extra arrow offset when creating buffers for TileDB write
    if (!arrow) {
        offsets.pop_back();
    }

    return {result, offsets};
}

template VarlenBufferPair to_varlen_buffers(
    std::vector<std::string>, bool arrow);

bool is_tiledb_uri(std::string_view uri) {
    return uri.find("tiledb://") == 0;
}

std::string rstrip_uri(std::string_view uri) {
    return std::regex_replace(std::string(uri), std::regex("/+$"), "");
}

std::optional<std::vector<uint8_t>> bitmap_to_uint8(
    const uint8_t* bitmap, size_t length, size_t offset) {
    if (bitmap == nullptr) {
        return std::nullopt;
    }

    std::vector<uint8_t> casted(length);
    ArrowBitsUnpackInt8(
        bitmap, offset, length, reinterpret_cast<int8_t*>(casted.data()));
    return casted;
}

std::shared_ptr<SOMAColumn> find_column_by_name(
    std::span<const std::shared_ptr<SOMAColumn>> columns,
    std::string_view name) {
    auto column_it = std::find_if(
        columns.begin(), columns.end(), [&](auto col) {
            return col->name() == name;
        });

    if (column_it == columns.end()) {
        throw TileDBSOMAError(fmt::format(
            "[ArrowAdapter][tiledb_schema_from_arrow_schema] Index column "
            "'{}' missing",
            name));
    }

    return *column_it;
}

std::string get_enmr_label(
    ArrowSchema* index_schema, ArrowSchema* value_schema) {
    std::string format(value_schema->format);
    format = (format == "u") ? "U" : (format == "z" ? "Z" : format);
    return std::string(index_schema->name) + "_" + format;
}

Enumeration get_enumeration(
    std::shared_ptr<Context> ctx,
    std::shared_ptr<Array> arr,
    ArrowSchema* index_schema,
    ArrowSchema* value_schema) {
    std::string new_way = util::get_enmr_label(index_schema, value_schema);
    std::string old_way = std::string(index_schema->name);
    try {
        // New-style names of the form {attr_name}_{arrow_format}, e.g. "foo_U"
        // for attributes written by tiledbsoma >= 1.16.0
        return ArrayExperimental::get_enumeration(*ctx, *arr, new_way);
    } catch (const std::exception& e) {
        // Old-style names of the form {attr_name}, e.g. "foo"
        // for attributes written by tiledbsoma < 1.16.0
        try {
            return ArrayExperimental::get_enumeration(*ctx, *arr, old_way);
        } catch (const std::exception& e) {
            throw TileDBSOMAError(fmt::format(
                "[get_enumeration] Could not find enumeration with name '{} or "
                "'{}'",
                new_way,
                old_way));
        }
    }
}

/**
 * Maps core Array/Group type enums to SOMA-style strings "SOMAArray" and
 * "SOMAGroup". Throws if the input value is neither one of those.
 */
std::string soma_type_from_tiledb_type(tiledb::Object::Type tiledb_type) {
    switch (tiledb_type) {
        case Object::Type::Array:
            return "SOMAArray";
        case Object::Type::Group:
            return "SOMAGroup";
        case Object::Type::Invalid:
            throw TileDBSOMAError(fmt::format(
                "[SOMAObject::open] Saw TileDB type Invalid ({}), which is "
                "neither Array ({}) nor Group ({})",
                static_cast<int>(tiledb_type),
                static_cast<int>(Object::Type::Array),
                static_cast<int>(Object::Type::Group)));
        default:
            throw TileDBSOMAError(fmt::format(
                "[SOMAObject::open] Saw unrecognized TileDB type ({}), which "
                "is neither Array ({}) nor Group ({}) nor Invalid ({})",
                static_cast<int>(tiledb_type),
                static_cast<int>(Object::Type::Array),
                static_cast<int>(Object::Type::Group),
                static_cast<int>(Object::Type::Invalid)));
    }
}

std::vector<std::array<DimType, 2>> get_fragment_non_empty_domain(
    const FragmentInfo& fragment_info, uint32_t idx, const Domain& domain) {
    std::vector<std::array<DimType, 2>> non_empy_domain(domain.ndim());

    for (size_t i = 0; i < domain.ndim(); ++i) {
        auto dimension = domain.dimension(i);

        switch (dimension.type()) {
            case TILEDB_UINT8: {
                uint8_t data[2];
                fragment_info.get_non_empty_domain(idx, i, data);
                non_empy_domain[i] = std::array<DimType, 2>({data[0], data[1]});
            } break;
            case TILEDB_INT8: {
                int8_t data[2];
                fragment_info.get_non_empty_domain(idx, i, data);
                non_empy_domain[i] = std::array<DimType, 2>({data[0], data[1]});
            } break;
            case TILEDB_UINT16: {
                uint16_t data[2];
                fragment_info.get_non_empty_domain(idx, i, data);
                non_empy_domain[i] = std::array<DimType, 2>({data[0], data[1]});
            } break;
            case TILEDB_INT16: {
                int16_t data[2];
                fragment_info.get_non_empty_domain(idx, i, data);
                non_empy_domain[i] = std::array<DimType, 2>({data[0], data[1]});
            } break;
            case TILEDB_UINT32: {
                uint32_t data[2];
                fragment_info.get_non_empty_domain(idx, i, data);
                non_empy_domain[i] = std::array<DimType, 2>({data[0], data[1]});
            } break;
            case TILEDB_INT32: {
                int32_t data[2];
                fragment_info.get_non_empty_domain(idx, i, data);
                non_empy_domain[i] = std::array<DimType, 2>({data[0], data[1]});
            } break;
            case TILEDB_UINT64: {
                uint64_t data[2];
                fragment_info.get_non_empty_domain(idx, i, data);
                non_empy_domain[i] = std::array<DimType, 2>({data[0], data[1]});
            } break;
            case TILEDB_INT64: {
                int64_t data[2];
                fragment_info.get_non_empty_domain(idx, i, data);
                non_empy_domain[i] = std::array<DimType, 2>({data[0], data[1]});
            } break;
            case TILEDB_FLOAT32: {
                float_t data[2];
                fragment_info.get_non_empty_domain(idx, i, data);
                non_empy_domain[i] = std::array<DimType, 2>({data[0], data[1]});
            } break;
            case TILEDB_FLOAT64: {
                double_t data[2];
                fragment_info.get_non_empty_domain(idx, i, data);
                non_empy_domain[i] = std::array<DimType, 2>({data[0], data[1]});
            } break;
            case TILEDB_STRING_ASCII: {
                auto data = fragment_info.non_empty_domain_var(idx, i);
                non_empy_domain[i] = std::array<DimType, 2>(
                    {data.first, data.second});
            } break;
            default:
                throw TileDBSOMAError(fmt::format(
                    "[get_fragment_non_empty_domain] Unsupported dimension "
                    "datatype {}",
                    tiledb::impl::type_to_str(dimension.type())));
        }
    }

    return non_empy_domain;
}
};  // namespace tiledbsoma::util
