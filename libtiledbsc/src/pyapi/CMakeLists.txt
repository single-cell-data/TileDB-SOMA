find_package(pybind11 REQUIRED)

pybind11_add_module(libtiledbsc
$<TARGET_OBJECTS:TILEDB_SC_OBJECTS>
libtiledbsc.cc
)

target_link_libraries(libtiledbsc
    PUBLIC
    tiledbsc
    TileDB::tiledb_shared
    )

target_include_directories(libtiledbsc
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src
    ${TILEDB_SC_EXPORT_HEADER_DIR}
)
