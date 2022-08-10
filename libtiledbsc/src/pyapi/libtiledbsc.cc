#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsc/tiledbsc>

#include "arrow_adapter.h"

#ifdef BUILD_COMMIT_HASH
#define VERSION BUILD_COMMIT_HASH
#else
#define VERSION "dev"
#endif

#define DENUM(x) .value(#x, TILEDB_##x)

using namespace tiledbsc;

namespace py = pybind11;
using namespace py::literals;

/**
 * @brief Convert ColumnBuffer to Arrow array.
 *
 * @param cb ColumnBuffer
 * @return py::object Arrow array
 */
py::object to_array(ColumnBuffer& cb) {
    auto pa = py::module::import("pyarrow");
    auto pa_array_import = pa.attr("Array").attr("_import_from_c");

    auto [array, schema] = ArrowAdapter::to_arrow(cb);
    return pa_array_import(py::capsule(array.get()), py::capsule(schema.get()));
}

/**
 * @brief Convert ArrayBuffers to Arrow table.
 *
 * @param cbs ArrayBuffers
 * @return py::object
 */
py::object to_table(ArrayBuffers& cbs) {
    auto pa = py::module::import("pyarrow");
    auto pa_table_from_arrays = pa.attr("Table").attr("from_arrays");

    py::list names;
    py::list arrays;

    for (auto& [name, column] : cbs) {
        names.append(name);
        arrays.append(to_array(*column));
    }

    return pa_table_from_arrays(arrays, names);
}

// TODO: Convert Arrow array to ColumnBuffer
/*
    _export_to_c(...) method of pyarrow.lib.StringArray instance
    Array._export_to_c(self, out_ptr, out_schema_ptr=0)

    Export to a C ArrowArray struct, given its pointer.

    If a C ArrowSchema struct pointer is also given, the array type
    is exported to it at the same time.

    Parameters
    ----------
    out_ptr: int
        The raw pointer to a C ArrowArray struct.
    out_schema_ptr: int (optional)
        The raw pointer to a C ArrowSchema struct.

    Be careful: if you don't pass the ArrowArray struct to a consumer,
    array memory will leak.  This is a low-level function intended for
    expert users.
*/

/**
 * @brief pybind11 bindings
 *
 */
PYBIND11_MODULE(libtiledbsc, m) {
    m.doc() = "SOMA acceleration library";

    m.def("version", []() {
        int major, minor, patch;
        tiledb_version(&major, &minor, &patch);
        return fmt::format(
            "libtiledbsc={} libtiledb={}.{}.{}", VERSION, major, minor, patch);
    });

    m.def(
        "config_logging",
        [](const std::string& level, const std::string& logfile) {
            LOG_CONFIG(level, logfile);
        },
        "level"_a,
        "logfile"_a = "");

    m.def("debug", &LOG_DEBUG, "message"_a = "");

    py::class_<SOMA>(m, "SOMA")
        .def(
            py::init(
                [](std::string_view uri,
                   std::optional<std::map<std::string, std::string>> config) {
                    if (config.has_value()) {
                        auto cfg = Config(*config);
                        return SOMA::open(uri, cfg);
                    } else {
                        return SOMA::open(uri);
                    }
                }),
            "uri"_a,
            "config"_a = py::none())
        .def("list_arrays", &SOMA::list_arrays)
        // SOMAQuery (0 = return value) will keep SOMA alive (1 = this)
        .def("query", &SOMA::query, py::keep_alive<0, 1>(), "name"_a = "soma");

    py::class_<SOMAQuery>(m, "SOMAQuery")
        .def("select_obs_attrs", &SOMAQuery::select_obs_attrs)
        .def("select_var_attrs", &SOMAQuery::select_var_attrs)
        .def("select_obs_ids", &SOMAQuery::select_obs_ids)
        .def("select_var_ids", &SOMAQuery::select_var_ids)
        .def("set_obs_condition", &SOMAQuery::set_obs_condition<uint64_t>)
        .def("set_obs_condition", &SOMAQuery::set_obs_condition<std::string>)
        .def("set_var_condition", &SOMAQuery::set_var_condition<uint64_t>)
        .def("set_var_condition", &SOMAQuery::set_var_condition<std::string>)
        .def(
            "next_results",
            [](SOMAQuery& sq)
                -> std::optional<std::map<std::string, py::object>> {
                auto buffers = sq.next_results();
                if (buffers.has_value()) {
                    std::map<std::string, py::object> results;
                    for (auto& [name, buffer] : *buffers) {
                        results[name] = to_table(buffer);
                    }
                    return results;
                }
                return std::nullopt;
            });

    py::class_<SOMACollection>(m, "SOMACollection")
        .def(
            py::init(
                [](std::string_view uri,
                   std::optional<std::map<std::string, std::string>> config) {
                    if (config.has_value()) {
                        auto cfg = Config(*config);
                        return SOMACollection::open(uri, cfg);
                    } else {
                        return SOMACollection::open(uri);
                    }
                }),
            "uri"_a,
            "config"_a = py::none())
        .def("list_somas", &SOMACollection::list_somas)
        // SOMACollection Query (0 = return value) will keep SOMACollection
        // alive (1 = this)
        .def("query", &SOMACollection::query, py::keep_alive<0, 1>());

    // clang-format off
    py::enum_<tiledb_query_condition_op_t>(m, "Condition", py::module_local()) 
        DENUM(LT)
        DENUM(LE)
        DENUM(GT)
        DENUM(GE)
        DENUM(EQ)
        DENUM(NE);
    // clang-format on
}
