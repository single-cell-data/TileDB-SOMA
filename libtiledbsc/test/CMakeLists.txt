############################################################
# Definitions
############################################################

# Set root source directory as compiler definition for use in tests
get_filename_component(TILEDBSC_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../" ABSOLUTE)
add_compile_definitions(TILEDBSC_SOURCE_ROOT="${TILEDBSC_SOURCE_ROOT}")

############################################################
# Dependencies
############################################################

find_package(TileDB_EP REQUIRED)

############################################################
# SC unit test
############################################################

if (TILEDBSC_TESTING)
    find_package(Catch_EP REQUIRED)

    add_executable(unit_sc EXCLUDE_FROM_ALL)
    target_link_libraries(unit_sc
        PUBLIC
        Catch2::Catch2WithMain
        TileDB::tiledb_shared
    )

    target_link_libraries(unit_sc
        PRIVATE
        $<TARGET_OBJECTS:buffer_set>
        $<TARGET_OBJECTS:managed_query>
    )

    target_include_directories(unit_sc
       PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/../src
        ../include)

    # Sources for tests
    target_sources(unit_sc PUBLIC
        unit_sc.cc
    )

    target_compile_definitions(unit_sc PRIVATE CATCH_CONFIG_MAIN)

    add_test(
        NAME "unit_sc"
        COMMAND $<TARGET_FILE:unit_sc>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_dependencies(
        unit_sc
        buffer_set
        managed_query
        tiledbsc
    )

    add_custom_target(build_tests)
    add_dependencies(build_tests
      unit_sc
      unit_util
      unit_buffer_set
      unit_query_result
      unit_ij_query
      unit_sc_arrow
      unit_managed_query
      # unit_py_query_result # Custom command, not a target
  )

    ############################################################
    # make check
    ############################################################

    # Add custom target 'check' in the root of the tree
    add_custom_target(
      check COMMAND ${CMAKE_CTEST_COMMAND} -V -C ${CMAKE_BUILD_TYPE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../
      DEPENDS
        build_tests
    )

endif()
