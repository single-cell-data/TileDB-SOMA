#include "common.h"

namespace tiledbsoma {

std::unordered_map<tiledb_datatype_t, std::string> _tdb_to_np_name_dtype = {
    {TILEDB_INT32, "int32"},
    {TILEDB_INT64, "int64"},
    {TILEDB_FLOAT32, "float32"},
    {TILEDB_FLOAT64, "float64"},
    {TILEDB_INT8, "int8"},
    {TILEDB_UINT8, "uint8"},
    {TILEDB_INT16, "int16"},
    {TILEDB_UINT16, "uint16"},
    {TILEDB_UINT32, "uint32"},
    {TILEDB_UINT64, "uint64"},
    {TILEDB_STRING_ASCII, "S"},
    {TILEDB_STRING_UTF8, "U1"},
    {TILEDB_CHAR, "S1"},
    {TILEDB_DATETIME_YEAR, "M8[Y]"},
    {TILEDB_DATETIME_MONTH, "M8[M]"},
    {TILEDB_DATETIME_WEEK, "M8[W]"},
    {TILEDB_DATETIME_DAY, "M8[D]"},
    {TILEDB_DATETIME_HR, "M8[h]"},
    {TILEDB_DATETIME_MIN, "M8[m]"},
    {TILEDB_DATETIME_SEC, "M8[s]"},
    {TILEDB_DATETIME_MS, "M8[ms]"},
    {TILEDB_DATETIME_US, "M8[us]"},
    {TILEDB_DATETIME_NS, "M8[ns]"},
    {TILEDB_DATETIME_PS, "M8[ps]"},
    {TILEDB_DATETIME_FS, "M8[fs]"},
    {TILEDB_DATETIME_AS, "M8[as]"},
    {TILEDB_TIME_HR, "m8[h]"},
    {TILEDB_TIME_MIN, "m8[m]"},
    {TILEDB_TIME_SEC, "m8[s]"},
    {TILEDB_TIME_MS, "m8[ms]"},
    {TILEDB_TIME_US, "m8[us]"},
    {TILEDB_TIME_NS, "m8[ns]"},
    {TILEDB_TIME_PS, "m8[ps]"},
    {TILEDB_TIME_FS, "m8[fs]"},
    {TILEDB_TIME_AS, "m8[as]"},
    {TILEDB_BLOB, "byte"},
    {TILEDB_BOOL, "bool"},
};

std::unordered_map<std::string, tiledb_datatype_t> _np_name_to_tdb_dtype = {
    {"int32", TILEDB_INT32},
    {"int64", TILEDB_INT64},
    {"float32", TILEDB_FLOAT32},
    {"float64", TILEDB_FLOAT64},
    {"int8", TILEDB_INT8},
    {"uint8", TILEDB_UINT8},
    {"int16", TILEDB_INT16},
    {"uint16", TILEDB_UINT16},
    {"uint32", TILEDB_UINT32},
    {"uint64", TILEDB_UINT64},
    {"datetime64[Y]", TILEDB_DATETIME_YEAR},
    {"datetime64[M]", TILEDB_DATETIME_MONTH},
    {"datetime64[W]", TILEDB_DATETIME_WEEK},
    {"datetime64[D]", TILEDB_DATETIME_DAY},
    {"datetime64[h]", TILEDB_DATETIME_HR},
    {"datetime64[m]", TILEDB_DATETIME_MIN},
    {"datetime64[s]", TILEDB_DATETIME_SEC},
    {"datetime64[ms]", TILEDB_DATETIME_MS},
    {"datetime64[us]", TILEDB_DATETIME_US},
    {"datetime64[ns]", TILEDB_DATETIME_NS},
    {"datetime64[ps]", TILEDB_DATETIME_PS},
    {"datetime64[fs]", TILEDB_DATETIME_FS},
    {"datetime64[as]", TILEDB_DATETIME_AS},
    /* duration types map to timedelta */
    {"timedelta64[h]", TILEDB_TIME_HR},
    {"timedelta64[m]", TILEDB_TIME_MIN},
    {"timedelta64[s]", TILEDB_TIME_SEC},
    {"timedelta64[ms]", TILEDB_TIME_MS},
    {"timedelta64[us]", TILEDB_TIME_US},
    {"timedelta64[ns]", TILEDB_TIME_NS},
    {"timedelta64[ps]", TILEDB_TIME_PS},
    {"timedelta64[fs]", TILEDB_TIME_FS},
    {"timedelta64[as]", TILEDB_TIME_AS},
    {"bool", TILEDB_BOOL},
};

py::dtype tdb_to_np_dtype(tiledb_datatype_t type, uint32_t cell_val_num) {
    if (type == TILEDB_CHAR || type == TILEDB_STRING_UTF8 ||
        type == TILEDB_STRING_ASCII) {
        std::string base_str = (type == TILEDB_STRING_UTF8) ? "|U" : "|S";
        if (cell_val_num < TILEDB_VAR_NUM)
            base_str += std::to_string(cell_val_num);
        return py::dtype(base_str);
    }

    if (cell_val_num == 1) {
        if (type == TILEDB_STRING_UTF16 || type == TILEDB_STRING_UTF32)
            TPY_ERROR_LOC("Unimplemented UTF16 or UTF32 string conversion!");
        if (type == TILEDB_STRING_UCS2 || type == TILEDB_STRING_UCS4)
            TPY_ERROR_LOC("Unimplemented UCS2 or UCS4 string conversion!");

        if (_tdb_to_np_name_dtype.count(type) == 1)
            return py::dtype(_tdb_to_np_name_dtype[type]);
    }

    if (cell_val_num == 2) {
        if (type == TILEDB_FLOAT32)
            return py::dtype("complex64");
        if (type == TILEDB_FLOAT64)
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
        "tiledb datatype not understood ('" + tiledb::impl::type_to_str(type) +
        "', cell_val_num: " + std::to_string(cell_val_num) + ")");
}

tiledb_datatype_t np_to_tdb_dtype(py::dtype type) {
    auto name = py::str(py::getattr(type, "name"));
    if (_np_name_to_tdb_dtype.count(name) == 1)
        return _np_name_to_tdb_dtype[name];

    auto kind = py::str(py::getattr(type, "kind"));
    if (kind.is(py::str("S")))
        return TILEDB_STRING_ASCII;
    if (kind.is(py::str("U")))
        return TILEDB_STRING_UTF8;

    TPY_ERROR_LOC("could not handle numpy dtype");
}

bool is_tdb_str(tiledb_datatype_t type) {
    switch (type) {
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_CHAR:
            return true;
        default:
            return false;
    }
}

/**
 * @brief Convert ArrayBuffers to Arrow table.
 *
 * @param cbs ArrayBuffers
 * @return py::object
 */
py::object _buffer_to_table(std::shared_ptr<ArrayBuffers> buffers) {
    auto pa = py::module::import("pyarrow");
    auto pa_table_from_arrays = pa.attr("Table").attr("from_arrays");
    auto pa_array_import = pa.attr("Array").attr("_import_from_c");
    auto pa_schema_import = pa.attr("Schema").attr("_import_from_c");

    py::list array_list;
    py::list names;

    for (auto& name : buffers->names()) {
        auto column = buffers->at(name);
        auto [pa_array, pa_schema] = ArrowAdapter::to_arrow(column);
        auto array = pa_array_import(
            py::capsule(pa_array.get()), py::capsule(pa_schema.get()));
        array_list.append(array);
        names.append(name);
    }

    return pa_table_from_arrays(array_list, names);
}

std::optional<py::object> to_table(
    std::optional<std::shared_ptr<ArrayBuffers>> buffers) {
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

        if (tdb_type == TILEDB_STRING_UTF8 || tdb_type == TILEDB_STRING_ASCII) {
            auto py_buf = py::array(py::dtype("|S1"), value_num, value);
            auto res = py_buf.attr("tobytes")().attr("decode")("UTF-8");
            results[py::str(key)] = res;
        } else {
            py::dtype value_type = tdb_to_np_dtype(tdb_type, 1);
            auto res = py::array(value_type, value_num, value).attr("item")(0);
            results[py::str(key)] = res;
        }
    }
    return results;
}

void set_metadata(
    SOMAObject& soma_object, const std::string& key, py::array value) {
    tiledb_datatype_t value_type = np_to_tdb_dtype(value.dtype());

    // For https://github.com/single-cell-data/TileDB-SOMA/pull/2900:
    // Ensure that all Python and R write paths use UTF-8 for string
    // metadata values.
    if (value_type == TILEDB_STRING_ASCII) {
        value_type = TILEDB_STRING_UTF8;
    }

    if (is_tdb_str(value_type) && value.size() > 1)
        throw py::type_error("array/list of strings not supported");

    py::buffer_info value_buffer = value.request();
    if (value_buffer.ndim != 1)
        throw py::type_error("Only 1D Numpy arrays can be stored as metadata");

    auto value_num = is_tdb_str(value_type) ? value.nbytes() : value.size();
    soma_object.set_metadata(
        key, value_type, value_num, value_num > 0 ? value.data() : nullptr);
}

}  // namespace tiledbsoma
