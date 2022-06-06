find_package(pybind11 REQUIRED)

##############################################################################
# module

pybind11_add_module(py_query_result_aux
  test/py_query_result_aux.cc
)

target_link_libraries(py_query_result_aux
    PUBLIC
    tiledbsc
    TileDB::tiledb_shared)

target_include_directories(py_query_result_aux
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    ${TILEDB_SC_EXPORT_HEADER_DIR}
)

##############################################################################
# find python

# Tell CMake to check the Python registry entry last on Windows
set(Python_FIND_REGISTRY "LAST")
# Tell CMake to prefer Python from the PATH
set(Python_FIND_STRATEGY "LOCATION")
find_package(Python COMPONENTS Interpreter Development REQUIRED)

##############################################################################
# test

add_custom_target(copy_py_test_file
  ALL DEPENDS
  ${CMAKE_CURRENT_BINARY_DIR}/test_py_query_result.py
)
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/test_py_query_result.py
  DEPENDS ${CMAKE_SOURCE_DIR}/src/lib/test/test_py_query_result.py
  COMMAND ${CMAKE_COMMAND} -E copy
          ${CMAKE_SOURCE_DIR}/src/lib/test/test_py_query_result.py
          ${CMAKE_CURRENT_BINARY_DIR}/test_py_query_result.py)

add_test(
    NAME "unit_py_query_result"
    COMMAND ${Python_EXECUTABLE} -m pytest test_py_query_result.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

#add_dependencies(unit_py_query_result copy_py_test_file)