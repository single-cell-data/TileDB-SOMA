#
# FindTileDB_EP.cmake
#
# Licensed under the MIT License.
# Copyright (c) TileDB, Inc.
#
# Finds the TileDB library, installing with an ExternalProject as necessary.

# If TileDB was installed as an EP, need to search the EP install path also.
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "${EP_INSTALL_PREFIX}")

if (FORCE_BUILD_TILEDB)
  find_package(TileDB CONFIG PATHS ${EP_INSTALL_PREFIX} NO_DEFAULT_PATH)
else()
  find_package(TileDB CONFIG)
endif()

if (TILEDB_FOUND)
  get_target_property(TILEDB_LIB TileDB::tiledb_shared IMPORTED_LOCATION_RELEASE)
  # If release build location not found, check for debug build location
  if (TILEDB_LIB MATCHES "NOTFOUND")
    get_target_property(TILEDB_LIB TileDB::tiledb_shared IMPORTED_LOCATION_DEBUG)
  endif()
  message(STATUS "Found TileDB: ${TILEDB_LIB}")
else()
  if (SUPERBUILD)
    message(STATUS "Adding TileDB as an external project")
    if (TILEDB_S3 STREQUAL "OFF")
      message(STATUS "TileDB will be built WITHOUT S3 support")
    endif()

    # Try to download prebuilt artifacts unless the user specifies to build from source.
    # NB When updating the pinned URLs here, please also update in file apis/r/tools/get_tarball.R
    if(DOWNLOAD_TILEDB_PREBUILT)
        if (WIN32) # Windows
          message(STATUS "Downloding TileDB windows-x86_64 binary")
          SET(DOWNLOAD_URL "https://github.com/TileDB-Inc/TileDB/releases/download/2.29.2/tiledb-windows-x86_64-2.29.2-2cd33d3.zip")
          SET(DOWNLOAD_SHA1 "5dbe87d66ea3840cb40aa96f7211a0e1b5c9d314")
        elseif(APPLE) # OSX
          if (CMAKE_OSX_ARCHITECTURES STREQUAL x86_64)
            message(STATUS "Downloding TileDB macos-x86_64 binary")
            SET(DOWNLOAD_URL "https://github.com/TileDB-Inc/TileDB/releases/download/2.29.2/tiledb-macos-x86_64-2.29.2-2cd33d3.tar.gz")
            SET(DOWNLOAD_SHA1 "ef68c7d86ae36b6126baf1ad824b4cae320a8ca1")
          elseif (CMAKE_OSX_ARCHITECTURES STREQUAL arm64)
            message(STATUS "Downloding TileDB macos-arm64 binary")
            SET(DOWNLOAD_URL "https://github.com/TileDB-Inc/TileDB/releases/download/2.29.2/tiledb-macos-arm64-2.29.2-2cd33d3.tar.gz")
            SET(DOWNLOAD_SHA1 "c9cdc34f91d0d9cc6374974c7f9ede5133ff820a")
          elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64)|(AMD64|amd64)|(^i.86$)")
            message(STATUS "Downloding TileDB macos-x86_64 binary")
            SET(DOWNLOAD_URL "https://github.com/TileDB-Inc/TileDB/releases/download/2.29.2/tiledb-macos-x86_64-2.29.2-2cd33d3.tar.gz")
            SET(DOWNLOAD_SHA1 "ef68c7d86ae36b6126baf1ad824b4cae320a8ca1")
          elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "^aarch64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "^arm")
            message(STATUS "Downloding TileDB macos-arm64 binary")
            SET(DOWNLOAD_URL "https://github.com/TileDB-Inc/TileDB/releases/download/2.29.2/tiledb-macos-arm64-2.29.2-2cd33d3.tar.gz")
            SET(DOWNLOAD_SHA1 "c9cdc34f91d0d9cc6374974c7f9ede5133ff820a")
          endif()

        else() # Linux
          if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64)")
            message(STATUS "Downloding TileDB linux-arm64 binary")
            SET(DOWNLOAD_URL "https://github.com/TileDB-Inc/TileDB/releases/download/2.29.2/tiledb-linux-arm64-2.29.2-2cd33d3.tar.gz")
            SET(DOWNLOAD_SHA1 "c175d860e0e4b42ae630f136724934c65c20bdd9")
          else()
            message(STATUS "Downloding TileDB linux-x86_64 binary")
            SET(DOWNLOAD_URL "https://github.com/TileDB-Inc/TileDB/releases/download/2.29.2/tiledb-linux-x86_64-2.29.2-2cd33d3.tar.gz")
            SET(DOWNLOAD_SHA1 "6c962b090c51506be9979294a21b65031a54e12e")
          endif()
        endif()

        ExternalProject_Add(ep_tiledb
          PREFIX "externals"
          URL ${DOWNLOAD_URL}
          URL_HASH SHA1=${DOWNLOAD_SHA1}
          CONFIGURE_COMMAND ""
          BUILD_COMMAND ""
          UPDATE_COMMAND ""
          PATCH_COMMAND ""
          TEST_COMMAND ""
          INSTALL_COMMAND
            ${CMAKE_COMMAND} -E copy_directory ${EP_BASE}/src/ep_tiledb ${EP_INSTALL_PREFIX}
          LOG_DOWNLOAD TRUE
          LOG_CONFIGURE FALSE
          LOG_BUILD FALSE
          LOG_INSTALL FALSE
        )
    else() # Build from source
        message(STATUS "Downloading TileDB source package")
        ExternalProject_Add(ep_tiledb
          PREFIX "externals"
          URL "https://github.com/TileDB-Inc/TileDB/archive/2.29.2.zip"
          URL_HASH SHA1=3b3652d40d4177522970f84f3b8f21d7eabf9962
          DOWNLOAD_NAME "tiledb.zip"
          CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${EP_INSTALL_PREFIX}
            -DCMAKE_PREFIX_PATH=${EP_INSTALL_PREFIX}
            -DTILEDB_S3=${TILEDB_S3}
            -DTILEDB_AZURE=${TILEDB_AZURE}
            -DTILEDB_GCS=${TILEDB_GCS}
            -DTILEDB_HDFS=${TILEDB_HDFS}
            -DTILEDB_SERIALIZATION=${TILEDB_SERIALIZATION}
            -DTILEDB_WERROR=${TILEDB_WERROR}
            -DTILEDB_REMOVE_DEPRECATIONS=${TILEDB_REMOVE_DEPRECATIONS}
            -DTILEDB_VERBOSE=${TILEDB_VERBOSE}
            -DTILEDB_TESTS=OFF
            -DCMAKE_BUILD_TYPE=$<CONFIG>
            -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
            -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
            -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            $<$<BOOL:${CMAKE_TOOLCHAIN_FILE}>:-DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}>
            -DCMAKE_POSITION_INDEPENDENT_CODE=ON
          UPDATE_COMMAND ""
          INSTALL_COMMAND
            ${CMAKE_COMMAND} --build . --target install
          LOG_DOWNLOAD TRUE
          LOG_CONFIGURE TRUE
          LOG_BUILD TRUE
          LOG_INSTALL TRUE
        )
    endif()

    list(APPEND FORWARD_EP_CMAKE_ARGS -DEP_TILEDB_BUILT=TRUE)
    list(APPEND EXTERNAL_PROJECTS ep_tiledb)
    message(STATUS "Finished adding TileDB external project")

  else()
    message(FATAL_ERROR "Unable to find TileDB library.")
  endif()
endif()

if (EP_TILEDB_BUILT AND TARGET TileDB::tiledb_shared)
  include(TileDBCommon)
  install_target_libs(TileDB::tiledb_shared)
endif()
