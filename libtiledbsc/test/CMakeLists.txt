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

    add_executable(unit_sc EXCLUDE_FROM_ALL
        $<TARGET_OBJECTS:TILEDB_SC_OBJECTS>
        unit_column_buffer.cc
        unit_managed_query.cc
    )

    target_link_libraries(unit_sc
      PRIVATE
        Catch2::Catch2WithMain
        TileDB::tiledb_shared
    )

    target_include_directories(unit_sc
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../src
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
    )

    target_compile_definitions(unit_sc PRIVATE CATCH_CONFIG_MAIN)

    add_test(
        NAME "unit_sc"
        COMMAND $<TARGET_FILE:unit_sc>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_dependencies(unit_sc
        tiledbsc
    )

    add_custom_target(build_tests)
    add_dependencies(build_tests
      unit_sc
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
