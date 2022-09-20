message(STATUS "")
message(STATUS "FILE ENTER libtiledbsoma/src/pyapi/CMakeLists.txt")
find_package(pybind11 REQUIRED)

pybind11_add_module(libtiledbsoma
    $<TARGET_OBJECTS:TILEDB_SOMA_OBJECTS>
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
    "XXX-INCN-1"
    ${CMAKE_CURRENT_SOURCE_DIR}
    "XXX-INCN-2"
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    "XXX-INCN-3"
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src
    "XXX-INCN-4"
    ${TILEDB_SOMA_EXPORT_HEADER_DIR}
    "XXX-INCN-5"
)
message(STATUS "FILE EXIT libtiledbsoma/src/pyapi/CMakeLists.txt")
message(STATUS "")
