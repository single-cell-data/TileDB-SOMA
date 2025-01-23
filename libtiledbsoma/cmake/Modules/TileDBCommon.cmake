#
# TileDBCommon.cmake
#
# Licensed under the MIT License.
# Copyright (c) TileDB, Inc.
#
# This file defines some common helper functions used by the external projects.
#

# Get library directory for multiarch linux distros
include(GNUInstallDirs)

#
# Stores the IMPORTED_LOCATION* target property of LIB_TARGET in RESULT_VAR.
# On Windows, preferentially tries the IMPORTED_IMPLIB* target property instead.
#
function(get_imported_location RESULT_VAR LIB_TARGET)
  if (WIN32)
    # Try several methods to find the imported location.
    get_target_property(TMP ${LIB_TARGET} IMPORTED_IMPLIB)
    if (TMP MATCHES "NOTFOUND")
      get_target_property(TMP ${LIB_TARGET} IMPORTED_IMPLIB_RELEASE)
    endif()
    if (TMP MATCHES "NOTFOUND")
      get_target_property(TMP ${LIB_TARGET} IMPORTED_IMPLIB_DEBUG)
    endif()
  endif()
  # Try several methods to find the imported location.
  if (TMP MATCHES "NOTFOUND" OR NOT WIN32)
    get_target_property(TMP ${LIB_TARGET} IMPORTED_LOCATION)
  endif()
  if (TMP MATCHES "NOTFOUND")
    get_target_property(TMP ${LIB_TARGET} IMPORTED_LOCATION_RELEASE)
  endif()
  if (TMP MATCHES "NOTFOUND")
    get_target_property(TMP ${LIB_TARGET} IMPORTED_LOCATION_DEBUG)
  endif()
  set(${RESULT_VAR} "${TMP}" PARENT_SCOPE)
endfunction()

#
# Concatenates the library name of the given imported target to the installation
# prefix and stores the result in RESULT_VAR. If the given library does not reside
# in the external projects build directory (i.e. it was not built by an EP), then
# just return the imported location.
#
function(get_installed_location RESULT_VAR LIB_TARGET)
  get_imported_location(IMP_LOC ${LIB_TARGET})
  if ("${IMP_LOC}" MATCHES "^${EP_BASE}")
    get_filename_component(LIB_NAME "${IMP_LOC}" NAME)
    set(INSTALLED_PATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${LIB_NAME}")
  else()
    set(INSTALLED_PATH "${IMP_LOC}")
  endif()
  set(${RESULT_VAR} "${INSTALLED_PATH}" PARENT_SCOPE)
endfunction()

#
# Adds imported libraries from the given target to the TileDB installation
# manifest.
#
function(install_target_libs LIB_TARGET)
  get_imported_location(TARGET_LIBRARIES ${LIB_TARGET})
  if (TARGET_LIBRARIES MATCHES "NOTFOUND")
    message(FATAL_ERROR "Could not determine library location for ${LIB_TARGET}")
  endif()

  get_filename_component(LIB_PATH ${TARGET_LIBRARIES} DIRECTORY)
  set(HEADERS_PATH ${LIB_PATH}/../include)

  if (WIN32 AND ${TARGET_LIBRARIES} MATCHES "${CMAKE_SHARED_LIBRARY_SUFFIX}$")
    install(FILES ${TARGET_LIBRARIES} DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(DIRECTORY ${HEADERS_PATH}/tiledb DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(FILES ${HEADERS_PATH}/tiledb_export.h DESTINATION ${CMAKE_INSTALL_BINDIR})
  else()
    get_filename_component(ABS_PATH ${TARGET_LIBRARIES} REALPATH)
    install(DIRECTORY ${LIB_PATH}/ DESTINATION lib)
    install(DIRECTORY ${HEADERS_PATH}/tiledb DESTINATION include)
    install(FILES ${HEADERS_PATH}/tiledb_export.h DESTINATION include)
  endif()
endfunction()
