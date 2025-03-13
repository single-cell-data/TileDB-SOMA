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
        // Note that this runs the release callbacks within the ArrowArray and
        // the ArrowSchema, freeing memory.
        auto array = pa_array_import(
            py::capsule(arrow_array->children[i]),
            py::capsule(arrow_schema->children[i]));
        array_list.append(array);

        // Already released: ensure there is no attempt at second free.
        arrow_array->children[i] = nullptr;
        arrow_schema->children[i] = nullptr;
    }
    // Already released: ensure there is no attempt at second free.
    arrow_array->n_children = 0;
    arrow_array->children = nullptr;
    arrow_schema->n_children = 0;
    arrow_schema->children = nullptr;

    return array_list;
}

void load_soma_array(py::module& m) {
    py::class_<SOMAArray, SOMAObject>(m, "SOMAArray")
        .def(
            py::init(
                [](std::string_view uri,
                   std::map<std::string, std::string> platform_config,
                   std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
                    return SOMAArray::open(
                        OpenMode::read,
                        uri,
                        std::make_shared<SOMAContext>(platform_config),
                        timestamp);
                }),
            "uri"_a,
            py::kw_only(),
            "platform_config"_a = py::dict(),
            "timestamp"_a = py::none())

        .def("__enter__", [](SOMAArray& array) { return array; })
        .def(
            "__exit__",
            [](SOMAArray& array,
               py::object exc_type,
               py::object exc_value,
               py::object traceback) { array.close(); })

        .def(
            "reopen",
            [](SOMAArray& array,
               OpenMode mode,
               std::optional<TimestampRange> timestamp) -> py::object {
                auto new_array = array.reopen(mode, timestamp);

                std::optional<std::string> soma_obj_type;
                try {
                    soma_obj_type = new_array->type();
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }

                if (!soma_obj_type) {
                    throw TileDBSOMAError(
                        "Unreachable code: The missing soma_object_type case "
                        "is already handled. This indicates an "
                        "unexpected failure to catch exceptions by "
                        "SOMAArray::reopen");
                }

                std::transform(
                    soma_obj_type->begin(),
                    soma_obj_type->end(),
                    soma_obj_type->begin(),
                    [](unsigned char c) { return std::tolower(c); });

                if (soma_obj_type == "somadataframe")
                    return py::cast(SOMADataFrame(*new_array));
                else if (soma_obj_type == "somapointclouddataframe")
                    return py::cast(SOMAPointCloudDataFrame(*new_array));
                else if (soma_obj_type == "somageometrydataframe")
                    return py::cast(SOMAGeometryDataFrame(*new_array));
                else if (soma_obj_type == "somasparsendarray")
                    return py::cast(SOMASparseNDArray(*new_array));
                else if (soma_obj_type == "somadensendarray")
                    return py::cast(SOMADenseNDArray(*new_array));

                throw TileDBSOMAError(
                    "Unreachable code: All possible SOMA object types are "
                    "already handled. This indicates a logic error or an "
                    "unexpected failure to catch exceptions by "
                    "SOMAArray::reopen");

                return py::none();  // Unreached, but appeases a compiler
                                    // warning
            },
            "mode"_a,
            "timestamp"_a = py::none())

        .def("close", &SOMAArray::close)
        .def_property_readonly(
            "closed",
            [](SOMAArray& array) -> bool { return not array.is_open(); })
        .def_property_readonly(
            "mode",
            [](SOMAArray& array) {
                return array.mode() == OpenMode::read ? "r" : "w";
            })
        .def_property_readonly(
            "schema",
            [](SOMAArray& array) -> py::object {
                auto pa = py::module::import("pyarrow");
                auto pa_schema_import = pa.attr("Schema").attr(
                    "_import_from_c");
                return pa_schema_import(
                    py::capsule(array.arrow_schema().get()));
            })
        .def("schema_config_options", &SOMAArray::schema_config_options)
        .def(
            "config_options_from_schema",
            &SOMAArray::config_options_from_schema)
        .def("context", &SOMAArray::ctx)

        .def("nnz", &SOMAArray::nnz, py::call_guard<py::gil_scoped_release>())

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
                ArrowTable arrow_table = array.get_non_empty_domain();
                return domainish_to_list(
                    arrow_table.first.get(), arrow_table.second.get());
            })

        .def(
            "domain",
            [](SOMAArray& array) {
                auto pa = py::module::import("pyarrow");
                ArrowTable arrow_table = array.get_soma_domain();
                return domainish_to_list(
                    arrow_table.first.get(), arrow_table.second.get());
            })

        .def(
            "maxdomain",
            [](SOMAArray& array) {
                auto pa = py::module::import("pyarrow");
                ArrowTable arrow_table = array.get_soma_maxdomain();
                return domainish_to_list(
                    arrow_table.first.get(), arrow_table.second.get());
            })

        .def(
            "non_empty_domain_slot",
            [](SOMAArray& array, std::string name, py::dtype dtype) {
                switch (np_to_tdb_dtype(dtype)) {
                    case TILEDB_UINT64:
                        return py::cast(
                            array.non_empty_domain_slot<uint64_t>(name));
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
                    case TILEDB_DATETIME_AS:
                    case TILEDB_INT64:
                        return py::cast(
                            array.non_empty_domain_slot<int64_t>(name));
                    case TILEDB_UINT32:
                        return py::cast(
                            array.non_empty_domain_slot<uint32_t>(name));
                    case TILEDB_INT32:
                        return py::cast(
                            array.non_empty_domain_slot<int32_t>(name));
                    case TILEDB_UINT16:
                        return py::cast(
                            array.non_empty_domain_slot<uint16_t>(name));
                    case TILEDB_INT16:
                        return py::cast(
                            array.non_empty_domain_slot<int16_t>(name));
                    case TILEDB_UINT8:
                        return py::cast(
                            array.non_empty_domain_slot<uint8_t>(name));
                    case TILEDB_INT8:
                        return py::cast(
                            array.non_empty_domain_slot<int8_t>(name));
                    case TILEDB_FLOAT64:
                        return py::cast(
                            array.non_empty_domain_slot<double>(name));
                    case TILEDB_FLOAT32:
                        return py::cast(
                            array.non_empty_domain_slot<float>(name));
                    case TILEDB_STRING_UTF8:
                    case TILEDB_STRING_ASCII:
                        return py::cast(
                            array.non_empty_domain_slot<std::string>(name));
                    default:
                        throw TileDBSOMAError(
                            "Unsupported dtype for nonempty domain.");
                }
            })

        .def(
            "non_empty_domain_slot_opt",
            [](SOMAArray& array, std::string name, py::dtype dtype) {
                switch (np_to_tdb_dtype(dtype)) {
                    case TILEDB_UINT64:
                        return py::cast(
                            array.non_empty_domain_slot_opt<uint64_t>(name));
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
                    case TILEDB_DATETIME_AS:
                    case TILEDB_INT64:
                        return py::cast(
                            array.non_empty_domain_slot_opt<int64_t>(name));
                    case TILEDB_UINT32:
                        return py::cast(
                            array.non_empty_domain_slot_opt<uint32_t>(name));
                    case TILEDB_INT32:
                        return py::cast(
                            array.non_empty_domain_slot_opt<int32_t>(name));
                    case TILEDB_UINT16:
                        return py::cast(
                            array.non_empty_domain_slot_opt<uint16_t>(name));
                    case TILEDB_INT16:
                        return py::cast(
                            array.non_empty_domain_slot_opt<int16_t>(name));
                    case TILEDB_UINT8:
                        return py::cast(
                            array.non_empty_domain_slot_opt<uint8_t>(name));
                    case TILEDB_INT8:
                        return py::cast(
                            array.non_empty_domain_slot_opt<int8_t>(name));
                    case TILEDB_FLOAT64:
                        return py::cast(
                            array.non_empty_domain_slot_opt<double>(name));
                    case TILEDB_FLOAT32:
                        return py::cast(
                            array.non_empty_domain_slot_opt<float>(name));
                    case TILEDB_STRING_UTF8:
                    case TILEDB_STRING_ASCII:
                        return py::cast(
                            array.non_empty_domain_slot_opt<std::string>(name));
                    default:
                        throw TileDBSOMAError(
                            "Unsupported dtype for nonempty domain.");
                }
            })

        .def(
            "soma_domain_slot",
            [](SOMAArray& array, std::string name, py::dtype dtype) {
                switch (np_to_tdb_dtype(dtype)) {
                    case TILEDB_UINT64:
                        return py::cast(array.soma_domain_slot<uint64_t>(name));
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
                    case TILEDB_DATETIME_AS:
                    case TILEDB_INT64:
                        return py::cast(array.soma_domain_slot<int64_t>(name));
                    case TILEDB_UINT32:
                        return py::cast(array.soma_domain_slot<uint32_t>(name));
                    case TILEDB_INT32:
                        return py::cast(array.soma_domain_slot<int32_t>(name));
                    case TILEDB_UINT16:
                        return py::cast(array.soma_domain_slot<uint16_t>(name));
                    case TILEDB_INT16:
                        return py::cast(array.soma_domain_slot<int16_t>(name));
                    case TILEDB_UINT8:
                        return py::cast(array.soma_domain_slot<uint8_t>(name));
                    case TILEDB_INT8:
                        return py::cast(array.soma_domain_slot<int8_t>(name));
                    case TILEDB_FLOAT64:
                        return py::cast(array.soma_domain_slot<double>(name));
                    case TILEDB_FLOAT32:
                        return py::cast(array.soma_domain_slot<float>(name));
                    case TILEDB_STRING_UTF8:
                    case TILEDB_STRING_ASCII: {
                        std::pair<std::string, std::string> str_domain;
                        return py::cast(std::make_pair("", ""));
                    }
                    default:
                        throw TileDBSOMAError(
                            "Unsupported dtype for Dimension's domain");
                }
            })

        .def(
            "soma_maxdomain_slot",
            [](SOMAArray& array, std::string name, py::dtype dtype) {
                switch (np_to_tdb_dtype(dtype)) {
                    case TILEDB_UINT64:
                        return py::cast(
                            array.soma_maxdomain_slot<uint64_t>(name));
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
                    case TILEDB_DATETIME_AS:
                    case TILEDB_INT64:
                        return py::cast(
                            array.soma_maxdomain_slot<int64_t>(name));
                    case TILEDB_UINT32:
                        return py::cast(
                            array.soma_maxdomain_slot<uint32_t>(name));
                    case TILEDB_INT32:
                        return py::cast(
                            array.soma_maxdomain_slot<int32_t>(name));
                    case TILEDB_UINT16:
                        return py::cast(
                            array.soma_maxdomain_slot<uint16_t>(name));
                    case TILEDB_INT16:
                        return py::cast(
                            array.soma_maxdomain_slot<int16_t>(name));
                    case TILEDB_UINT8:
                        return py::cast(
                            array.soma_maxdomain_slot<uint8_t>(name));
                    case TILEDB_INT8:
                        return py::cast(
                            array.soma_maxdomain_slot<int8_t>(name));
                    case TILEDB_FLOAT64:
                        return py::cast(
                            array.soma_maxdomain_slot<double>(name));
                    case TILEDB_FLOAT32:
                        return py::cast(array.soma_maxdomain_slot<float>(name));
                    case TILEDB_STRING_UTF8:
                    case TILEDB_STRING_ASCII: {
                        std::pair<std::string, std::string> str_domain;
                        return py::cast(std::make_pair("", ""));
                    }
                    default:
                        throw TileDBSOMAError(
                            "Unsupported dtype for Dimension's domain");
                }
            })

        .def_property_readonly("dimension_names", &SOMAArray::dimension_names)

        .def(
            "consolidate_and_vacuum",
            &SOMAArray::consolidate_and_vacuum,
            py::arg(
                "modes") = std::vector<std::string>{"fragment_meta", "commits"})

        .def_property_readonly(
            "meta",
            [](SOMAArray& array) -> py::dict {
                return meta(array.get_metadata());
            })

        .def(
            "set_metadata",
            set_metadata,
            py::arg("key"),
            py::arg("value"),
            py::arg("force") = false)

        .def(
            "delete_metadata",
            &SOMAArray::delete_metadata,
            py::arg("key"),
            py::arg("force") = false)

        .def(
            "get_metadata",
            py::overload_cast<const std::string&>(&SOMAArray::get_metadata))

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
        .def_property_readonly(
            "tiledbsoma_has_upgraded_shape", &SOMAArray::has_current_domain)

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
                    return array.can_upgrade_shape(
                        newshape, "tiledbsoma_can_upgrade_shape");
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
