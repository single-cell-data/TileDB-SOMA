/**
 * @file  common.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines common functions for the SOMA PyBind layer.
 */

#include "common.h"

using namespace pybind11::literals;  // to bring in the `_a` literal

namespace tiledbsoma {

DomainRange encode_domain(std::string_view format, py::object domain) {
    auto encode_element = [&]<typename T>() -> DomainRange {
        if (domain.is_none()) {
            return std::optional<std::pair<T, T>>();
        }
        return std::make_optional<std::pair<T, T>>(domain.cast<std::pair<T, T>>());
    };

    switch (tiledbsoma::common::arrow::to_tiledb_format(format)) {
        case TILEDB_UINT8:
            return encode_element.template operator()<uint8_t>();
        case TILEDB_UINT16:
            return encode_element.template operator()<uint16_t>();
        case TILEDB_UINT32:
            return encode_element.template operator()<uint32_t>();
        case TILEDB_UINT64:
            return encode_element.template operator()<uint64_t>();
        case TILEDB_INT8:
            return encode_element.template operator()<int8_t>();
        case TILEDB_INT16:
            return encode_element.template operator()<int16_t>();
        case TILEDB_INT32:
            return encode_element.template operator()<int32_t>();
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_INT64:
            return encode_element.template operator()<int64_t>();
        case TILEDB_FLOAT32:
            return encode_element.template operator()<float_t>();
        case TILEDB_FLOAT64:
            return encode_element.template operator()<double_t>();
        case TILEDB_CHAR:
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_BLOB:
        case TILEDB_GEOM_WKT:
        case TILEDB_GEOM_WKB:
            return encode_element.template operator()<std::string>();
        default:
            throw py::type_error(
                "[encode_domain] Unsupported type " +
                tiledb::impl::type_to_str(tiledbsoma::common::arrow::to_tiledb_format(format)));
    }
}

std::unordered_map<common::DataType, std::string> _tdb_to_np_name_dtype = {
    {common::DataType::INT32, "int32"},
    {common::DataType::INT64, "int64"},
    {common::DataType::FLOAT32, "float32"},
    {common::DataType::FLOAT64, "float64"},
    {common::DataType::INT8, "int8"},
    {common::DataType::UINT8, "uint8"},
    {common::DataType::INT16, "int16"},
    {common::DataType::UINT16, "uint16"},
    {common::DataType::UINT32, "uint32"},
    {common::DataType::UINT64, "uint64"},
    {common::DataType::STRING_ASCII, "S"},
    {common::DataType::STRING_UTF8, "U1"},
    {common::DataType::CHAR, "S1"},
    {common::DataType::DATETIME_YEAR, "M8[Y]"},
    {common::DataType::DATETIME_MONTH, "M8[M]"},
    {common::DataType::DATETIME_WEEK, "M8[W]"},
    {common::DataType::DATETIME_DAY, "M8[D]"},
    {common::DataType::DATETIME_HR, "M8[h]"},
    {common::DataType::DATETIME_MIN, "M8[m]"},
    {common::DataType::DATETIME_SEC, "M8[s]"},
    {common::DataType::DATETIME_MS, "M8[ms]"},
    {common::DataType::DATETIME_US, "M8[us]"},
    {common::DataType::DATETIME_NS, "M8[ns]"},
    {common::DataType::DATETIME_PS, "M8[ps]"},
    {common::DataType::DATETIME_FS, "M8[fs]"},
    {common::DataType::DATETIME_AS, "M8[as]"},
    {common::DataType::TIME_HR, "m8[h]"},
    {common::DataType::TIME_MIN, "m8[m]"},
    {common::DataType::TIME_SEC, "m8[s]"},
    {common::DataType::TIME_MS, "m8[ms]"},
    {common::DataType::TIME_US, "m8[us]"},
    {common::DataType::TIME_NS, "m8[ns]"},
    {common::DataType::TIME_PS, "m8[ps]"},
    {common::DataType::TIME_FS, "m8[fs]"},
    {common::DataType::TIME_AS, "m8[as]"},
    {common::DataType::BLOB, "V"},
    {common::DataType::BOOL, "bool"},
};

std::unordered_map<std::string, common::DataType> _np_name_to_tdb_dtype = {
    {"int32", common::DataType::INT32},
    {"int64", common::DataType::INT64},
    {"float32", common::DataType::FLOAT32},
    {"float64", common::DataType::FLOAT64},
    {"int8", common::DataType::INT8},
    {"uint8", common::DataType::UINT8},
    {"int16", common::DataType::INT16},
    {"uint16", common::DataType::UINT16},
    {"uint32", common::DataType::UINT32},
    {"uint64", common::DataType::UINT64},
    {"datetime64[Y]", common::DataType::DATETIME_YEAR},
    {"datetime64[M]", common::DataType::DATETIME_MONTH},
    {"datetime64[W]", common::DataType::DATETIME_WEEK},
    {"datetime64[D]", common::DataType::DATETIME_DAY},
    {"datetime64[h]", common::DataType::DATETIME_HR},
    {"datetime64[m]", common::DataType::DATETIME_MIN},
    {"datetime64[s]", common::DataType::DATETIME_SEC},
    {"datetime64[ms]", common::DataType::DATETIME_MS},
    {"datetime64[us]", common::DataType::DATETIME_US},
    {"datetime64[ns]", common::DataType::DATETIME_NS},
    {"datetime64[ps]", common::DataType::DATETIME_PS},
    {"datetime64[fs]", common::DataType::DATETIME_FS},
    {"datetime64[as]", common::DataType::DATETIME_AS},
    /* duration types map to timedelta */
    {"timedelta64[h]", common::DataType::TIME_HR},
    {"timedelta64[m]", common::DataType::TIME_MIN},
    {"timedelta64[s]", common::DataType::TIME_SEC},
    {"timedelta64[ms]", common::DataType::TIME_MS},
    {"timedelta64[us]", common::DataType::TIME_US},
    {"timedelta64[ns]", common::DataType::TIME_NS},
    {"timedelta64[ps]", common::DataType::TIME_PS},
    {"timedelta64[fs]", common::DataType::TIME_FS},
    {"timedelta64[as]", common::DataType::TIME_AS},
    {"bool", common::DataType::BOOL},
};

py::dtype tdb_to_np_dtype(common::DataType type, uint32_t cell_val_num) {
    if (type == common::DataType::BLOB) {
        std::string base_str = "|V";
        if (cell_val_num < TILEDB_VAR_NUM)
            base_str += std::to_string(cell_val_num);
        return py::dtype(base_str);
    }

    if (type == common::DataType::CHAR || type == common::DataType::STRING_UTF8 ||
        type == common::DataType::STRING_ASCII) {
        std::string base_str = (type == common::DataType::STRING_UTF8) ? "|U" : "|S";
        if (cell_val_num < TILEDB_VAR_NUM)
            base_str += std::to_string(cell_val_num);
        return py::dtype(base_str);
    }

    if (cell_val_num == 1) {
        if (type == common::DataType::STRING_UTF16 || type == common::DataType::STRING_UTF32)
            TPY_ERROR_LOC("Unimplemented UTF16 or UTF32 string conversion!");
        if (type == common::DataType::STRING_UCS2 || type == common::DataType::STRING_UCS4)
            TPY_ERROR_LOC("Unimplemented UCS2 or UCS4 string conversion!");

        if (_tdb_to_np_name_dtype.count(type) == 1)
            return py::dtype(_tdb_to_np_name_dtype[type]);
    }

    if (cell_val_num == 2) {
        if (type == common::DataType::FLOAT32)
            return py::dtype("complex64");
        if (type == common::DataType::FLOAT64)
            return py::dtype("complex128");
    }

    if (cell_val_num == TILEDB_VAR_NUM)
        return tdb_to_np_dtype(type, 1);

    if (cell_val_num > 1) {
        py::dtype base_dtype = tdb_to_np_dtype(type, 1);
        py::tuple rec_elem = py::make_tuple("", base_dtype);
        py::list rec_list;
        for (size_t i = 0; i < cell_val_num; i++)
            rec_list.append(rec_elem);
        // note: we call the 'dtype' constructor b/c py::dtype does not accept
        // list
        auto np = py::module::import("numpy");
        auto np_dtype = np.attr("dtype");
        return np_dtype(rec_list);
    }

    TPY_ERROR_LOC(
        "tiledb datatype not understood ('" + std::string(common::getName(type)) +
        "', cell_val_num: " + std::to_string(cell_val_num) + ")");
}

common::DataType np_to_tdb_dtype(py::dtype type) {
    auto name = py::str(py::getattr(type, "name"));
    if (_np_name_to_tdb_dtype.count(name) == 1)
        return _np_name_to_tdb_dtype[name];

    auto kind = py::str(py::getattr(type, "kind"));

    if (kind.is(py::str("S")))
        return common::DataType::STRING_UTF8;
    // Numpy encodes strings as UTF-32
    if (kind.is(py::str("U")))
        TPY_ERROR_LOC("[np_to_tdb_dtype] UTF-32 encoded strings are not supported");

    // No std::format in C++17, and, including logger/fmt headers
    // is tetchy here.
    std::stringstream ss;
    ss << "[np_to_tdb_dtype] Could not handle numpy dtype of kind '";
    ss << kind.operator std::string();
    ss << "'";
    TPY_ERROR_LOC(ss.str());
}

bool is_tdb_str(common::DataType type) {
    switch (type) {
        case common::DataType::STRING_ASCII:
        case common::DataType::STRING_UTF8:
        case common::DataType::CHAR:
            return true;
        default:
            return false;
    }
}

/**
 * @brief Convert ArrayBuffers to Arrow table.
 *
 * @param buffers ArrayBuffers
 * @return py::object pa.Table
 */
py::object _buffer_to_table(std::shared_ptr<common::ArrayBuffers> buffers) {
    auto pa = py::module::import("pyarrow");
    auto pa_array_import = pa.attr("Array").attr("_import_from_c");
    auto pa_dtype_import = pa.attr("DataType").attr("_import_from_c");
    auto pa_table_from_arrays = pa.attr("Table").attr("from_arrays");

    py::list array_list;
    py::list field_list;

    auto column_names = buffers->names();

    auto arrays = ArrowAdapter::buffer_to_arrow(buffers);

    for (size_t i = 0; i < column_names.size(); ++i) {
        auto& [pa_array, pa_schema] = arrays[i];
        auto nullable = (pa_schema->flags & ARROW_FLAG_NULLABLE) != 0;
        auto dtype = pa_dtype_import(py::capsule(pa_schema.get()));
        array_list.append(pa_array_import(py::capsule(pa_array.get()), dtype));
        field_list.append(pa.attr("field")(column_names[i], dtype, nullable));
    }

    return pa_table_from_arrays(array_list, "schema"_a = pa.attr("schema")(field_list));
}

std::optional<py::object> to_table(std::optional<std::shared_ptr<common::ArrayBuffers>> buffers) {
    // If more data was read, convert it to an arrow table and return
    if (buffers.has_value()) {
        return _buffer_to_table(*buffers);
    }

    // No data was read, the query is complete, return nullopt
    return std::nullopt;
}

py::dict meta(std::map<std::string, MetadataValue> metadata_mapping) {
    py::dict results;

    for (auto [key, val] : metadata_mapping) {
        auto [tdb_type, value_num, value] = val;

        if (tdb_type == common::DataType::STRING_UTF8 || tdb_type == common::DataType::STRING_ASCII) {
            // Empty strings stored as nullptr have a value_num of 1 and a \x00
            // value
            if (value_num == 1 && value == nullptr) {
                results[py::str(key)] = "";
            } else {
                auto py_buf = py::array(py::dtype("|S1"), value_num, value);
                results[py::str(key)] = py_buf.attr("tobytes")().attr("decode")("UTF-8");
            }
        } else if (tdb_type == common::DataType::BLOB) {
            py::dtype value_type = tdb_to_np_dtype(tdb_type, value_num);
            results[py::str(key)] = py::array(value_type, value_num, value).attr("item")(0);
        } else {
            py::dtype value_type = tdb_to_np_dtype(tdb_type, value_num);
            results[py::str(key)] = py::array(value_type, value_num, value).attr("item")(0);
        }
    }
    return results;
}

void set_metadata(SOMAObject& soma_object, const std::string& key, py::array value, bool force) {
    common::DataType value_type = np_to_tdb_dtype(value.dtype());

    // For https://github.com/single-cell-data/TileDB-SOMA/pull/2900:
    // Ensure that all Python and R write paths use UTF-8 for string
    // metadata values.
    if (value_type == common::DataType::STRING_ASCII) {
        value_type = common::DataType::STRING_UTF8;
    }

    if (is_tdb_str(value_type) && value.size() > 1)
        throw py::type_error("array/list of strings not supported");

    py::buffer_info value_buffer = value.request();
    if (value_buffer.ndim != 1)
        throw py::type_error("Only 1D Numpy arrays can be stored as metadata");

    auto value_num = is_tdb_str(value_type) ? value.nbytes() : value.size();

    if (is_tdb_str(value_type) && value_num > 0) {
        // If an empty string is passed by default results in a NULL byte
        switch (value_type) {
            case common::DataType::STRING_UTF8:
                value_num = sanitize_string(
                    std::span<const uint8_t>(static_cast<const uint8_t*>(value.data()), value_num), value_num);

                break;
            default:
                // No std::format in C++17, and, including logger/fmt headers
                // is tetchy here.
                std::stringstream ss;
                ss << "[set_metadata] Unsupported string encoding '" << common::getName(value_type) << "' for key '"
                   << key << "'";
                throw TileDBSOMAError(ss.str());
        }
    }

    // `sanitize_string` will return 0 if the firts byte is NULL and the string
    // is length 1 to account for empty numpy ascii strings. But for keys this
    // will never be the case since empty python strings are correctly converter
    // to std::string. So detecting a difference in the passed key length and
    // the sanitized length is an error.

    // Moreover all python strings are utf-8 encoded so eny NULL byte included
    // in the bytestream will be indeed a NULL byte
    if (sanitize_string(std::span<const char>(key.c_str(), key.length()), key.length()) != key.length()) {
        throw TileDBSOMAError("[set_metadata] Key contains NULL bytes");
    }

    soma_object.set_metadata(
        key,
        value_type,
        value_num,
        value_num > 0 ? value.data() : nullptr,
        force);  // The force flag is intended to only be toggled for testing in
                 // fault-injection scenarios
}

}  // namespace tiledbsoma
