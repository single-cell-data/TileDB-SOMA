find_package(pybind11 REQUIRED)

pybind11_add_module(libtiledbsoma
    $<TARGET_OBJECTS:TILEDBSOMA_OBJECTS>
    libtiledbsoma.cc
)

target_compile_definitions(libtiledbsoma PRIVATE
    -DBUILD_COMMIT_HASH="${BUILD_COMMIT_HASH}"
)

target_link_libraries(libtiledbsoma
    PUBLIC
    TileDB::tiledb_shared
)

target_include_directories(libtiledbsoma
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src
    ${TILEDBSOMA_EXPORT_HEADER_DIR}
)
