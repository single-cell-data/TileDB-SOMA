find_package(pybind11 REQUIRED)

pybind11_add_module(pytiledbsc
  pytiledbsc.cc
)

target_link_libraries(pytiledbsc
    PUBLIC
    tiledbsc
    TileDB::tiledb_shared
    )

target_include_directories(pytiledbsc
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    ${TILEDB_SC_EXPORT_HEADER_DIR}
)