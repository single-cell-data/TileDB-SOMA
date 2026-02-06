/**
 * @file   soma_array.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the SOMAArray bindings.
 */

#include "common.h"

#define DENUM(x) .value(#x, TILEDB_##x)
namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

// This is shared code for non-empty domain, domain, and maxdomain.  Returns a
// Python list of Arrow arrays. Python code can convert these further.
py::list domainish_to_list(ArrowArray* arrow_array, ArrowSchema* arrow_schema) {
    auto pa = py::module::import("pyarrow");
    auto pa_array_import = pa.attr("Array").attr("_import_from_c");

    py::list array_list;
    for (int i = 0; i < arrow_array->n_children; i++) {
        // "_import_from_c" implements Arrow move semantics, meaning
        // release is set to NULL for each array imported. Release
        // should be called on the parent array.
        auto array = pa_array_import(py::capsule(arrow_array->children[i]), py::capsule(arrow_schema->children[i]));
        array_list.append(array);
    }

    return array_list;
}

void load_soma_array(py::module& m) {
    py::class_<SOMAArray, SOMAObject>(m, "SOMAArray")
        .def(
            py::init([](std::string_view uri,
                        std::map<std::string, std::string> platform_config,
                        std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
                return SOMAArray::open(
                    OpenMode::soma_read, uri, std::make_shared<SOMAContext>(platform_config), timestamp);
            }),
            "uri"_a,
            py::kw_only(),
            "platform_config"_a = py::dict(),
            "timestamp"_a = py::none())

        .def("__enter__", [](SOMAArray& array) { return array; })
        .def(
            "__exit__",
            [](SOMAArray& array, py::object exc_type, py::object exc_value, py::object traceback) { array.close(); })

        .def_property_readonly("type", &SOMAArray::type)
        .def("close", &SOMAArray::close)
        .def_property_readonly("closed", [](SOMAArray& array) -> bool { return not array.is_open(); })
        .def_property_readonly(
            "mode",
            [](SOMAArray& array) {
                OpenMode soma_mode = array.mode();
                switch (soma_mode) {
                    case OpenMode::soma_read:
                        return "r";
                    case OpenMode::soma_write:
                        return "w";
                    case OpenMode::soma_delete:
                        return "d";
                    default:
                        throw TileDBSOMAError("Internal error: unrecognized mode.");
                }
            })
        .def_property_readonly(
            "schema",
            [](SOMAArray& array) -> py::object {
                auto pa = py::module::import("pyarrow");
                auto pa_schema_import = pa.attr("Schema").attr("_import_from_c");

                try {
                    return pa_schema_import(py::capsule(array.arrow_schema().get()));
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
            })
        .def("schema_config_options", &SOMAArray::schema_config_options)
        .def("context", &SOMAArray::ctx)

        .def(
            "nnz",
            [](SOMAArray& array) {
                try {
                    py::gil_scoped_release release;
                    auto retval = array.nnz();
                    py::gil_scoped_acquire acquire;
                    return retval;
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
            })

        .def(
            "fragment_cell_count",
            [](SOMAArray& array) {
                try {
                    py::gil_scoped_release release;
                    auto retval = array.fragment_cell_count();
                    py::gil_scoped_acquire acquire;
                    return retval;
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
            })

        .def_property_readonly("uri", &SOMAArray::uri)

        .def_property_readonly(
            "timestamp",
            [](SOMAArray& array) -> py::object {
                if (!array.timestamp().has_value())
                    return py::none();
                return py::cast(array.timestamp()->second);
            })

        .def("domainish_to_list", domainish_to_list)

        .def(
            "non_empty_domain",
            [](SOMAArray& array) {
                common::arrow::ArrowTable arrow_table = array.get_non_empty_domain();
                return domainish_to_list(arrow_table.first.get(), arrow_table.second.get());
            })

        .def(
            "domain",
            [](SOMAArray& array) {
                auto pa = py::module::import("pyarrow");
                common::arrow::ArrowTable arrow_table = array.get_soma_domain();
                return domainish_to_list(arrow_table.first.get(), arrow_table.second.get());
            })

        .def(
            "maxdomain",
            [](SOMAArray& array) {
                auto pa = py::module::import("pyarrow");
                common::arrow::ArrowTable arrow_table = array.get_soma_maxdomain();
                return domainish_to_list(arrow_table.first.get(), arrow_table.second.get());
            })

        .def_property_readonly("dimension_names", &SOMAArray::dimension_names)

        .def(
            "consolidate_and_vacuum",
            &SOMAArray::consolidate_and_vacuum,
            py::arg("modes") = std::vector<std::string>{"fragment_meta", "commits"})

        .def_property_readonly("meta", [](SOMAArray& array) -> py::dict { return meta(array.get_metadata()); })

        .def("set_metadata", set_metadata, py::arg("key"), py::arg("value"), py::arg("force") = false)

        .def("delete_metadata", &SOMAArray::delete_metadata, py::arg("key"), py::arg("force") = false)

        .def("get_metadata", py::overload_cast<const std::string&>(&SOMAArray::get_metadata))

        .def("has_metadata", &SOMAArray::has_metadata)

        .def("metadata_num", &SOMAArray::metadata_num)

        // These are for SparseNDArray and DenseNDArray both:
        // * tiledbsoma_has_upgraded_shape
        // * resize
        // * can_resize
        // * tiledbsoma_upgrade_shape
        // * tiledbsoma_can_upgrade_shape
        // We don't have CommonNDArray base class in pybind11, and it's probably
        // not worth it.  These are exposed to the user-facing API only for
        // SparseNDArray and DenseNDArray and not for DataFrame.
        .def_property_readonly("tiledbsoma_has_upgraded_shape", &SOMAArray::has_current_domain)

        .def(
            "resize",
            [](SOMAArray& array, const std::vector<int64_t>& newshape) {
                try {
                    array.resize(newshape, "resize");
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "newshape"_a)
        .def(
            "can_resize",
            [](SOMAArray& array, const std::vector<int64_t>& newshape) {
                try {
                    return array.can_resize(newshape, "can_resize");
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "newshape"_a)

        .def(
            "tiledbsoma_upgrade_shape",
            [](SOMAArray& array, const std::vector<int64_t>& newshape) {
                try {
                    array.upgrade_shape(newshape, "tiledbsoma_upgrade_shape");
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "newshape"_a)
        .def(
            "tiledbsoma_can_upgrade_shape",
            [](SOMAArray& array, const std::vector<int64_t>& newshape) {
                try {
                    return array.can_upgrade_shape(newshape, "tiledbsoma_can_upgrade_shape");
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "newshape"_a)
        .def(
            "get_column",
            [](SOMAArray& array, const std::string name) {
                try {
                    return array.get_column(name);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "name"_a);
}
}  // namespace libtiledbsomacpp
