if(DEFINED ENV{LIBTILEDBSOMA_DEBUG_BUILD})
  message(STATUS "")
  message(STATUS "START libtiledbsoma/src/pyapi/CMakeLists.txt")
endif()
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
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src
    ${TILEDB_SOMA_EXPORT_HEADER_DIR}
)
if(DEFINED ENV{LIBTILEDBSOMA_DEBUG_BUILD})
  message(STATUS "END libtiledbsoma/src/pyapi/CMakeLists.txt")
  message(STATUS "")
endif()
