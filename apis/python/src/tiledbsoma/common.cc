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
      throw TileDBSOMAError("Unimplemented UTF16 or UTF32 string conversion!");
    if (type == TILEDB_STRING_UCS2 || type == TILEDB_STRING_UCS4)
      throw TileDBSOMAError("Unimplemented UCS2 or UCS4 string conversion!");

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

  throw TileDBSOMAError("tiledb datatype not understood ('" +
        tiledb::impl::type_to_str(type) +
        "', cell_val_num: " + std::to_string(cell_val_num) + ")");
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
        auto [pa_array, pa_schema] = ArrowAdapter::to_arrow(column, true);
        auto array = pa_array_import(py::capsule(pa_array.get()), 
                                     py::capsule(pa_schema.get()));
        array_list.append(array);
        names.append(name);
    }

    return pa_table_from_arrays(array_list, names);
}

std::optional<py::object> to_table(
    std::optional<std::shared_ptr<ArrayBuffers>> buffers){
    // If more data was read, convert it to an arrow table and return
    if (buffers.has_value()) {
        return _buffer_to_table(*buffers);
    }

    // No data was read, the query is complete, return nullopt
    return std::nullopt;
}

}