# https://pybind11.readthedocs.io/en/stable/compiling.html#findpython-mode
#find_package(Python COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(libtiledbsoma
    $<TARGET_OBJECTS:TILEDB_SOMA_OBJECTS>
    libtiledbsoma.cc
)

# Update list in parent scope
list(APPEND TILEDBSOMA_INSTALL_TARGETS libtiledbsoma)
set(TILEDBSOMA_INSTALL_TARGETS ${TILEDBSOMA_INSTALL_TARGETS} PARENT_SCOPE)

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
    ${CMAKE_CURRENT_SOURCE_DIR}/../../external/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src
    ${TILEDB_SOMA_EXPORT_HEADER_DIR}
    $<TARGET_PROPERTY:spdlog::spdlog,INTERFACE_INCLUDE_DIRECTORIES>
)

# Set rpath to be relative to the .so.
# This allows tiledbsoma shared objects to find the tiledb shared object in the
# same directory, for example in the installed python module directory.
if (APPLE)
  set_target_properties(libtiledbsoma PROPERTIES INSTALL_RPATH "@loader_path/")
else()
  set_target_properties(libtiledbsoma PROPERTIES INSTALL_RPATH "$ORIGIN/")
endif()
