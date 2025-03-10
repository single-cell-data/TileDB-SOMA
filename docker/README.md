# Building/Installing TileDB-SOMA in Docker images

Sample Dockerfiles that `pip install tiledbsoma` (including installing required system deps).

## Debugging Tips

### `pip install --no-clean tiledbsoma` <a id="no-clean"></a>
This keeps around various logs (example below), when a `pip install` fails.

### [vcpkg] requires `VCPKG_FORCE_SYSTEM_BINARIES=1` on ARM <a id="vcpkg"></a>

Building on ARM without `VCPKG_FORCE_SYSTEM_BINARIES=1` fails like:

```bash
make install
# Building Release build
# -- Using default install prefix /TileDB-SOMA/dist. To control CMAKE_INSTALL_PREFIX, set OVERRIDE_INSTALL_PREFIX=OFF
# -- Install prefix is /TileDB-SOMA/dist.
# -- The C compiler identification is GNU 13.3.0
# -- The CXX compiler identification is GNU 13.3.0
# -- Detecting C compiler ABI info
# -- Detecting C compiler ABI info - done
# -- Check for working C compiler: /usr/bin/cc - skipped
# -- Detecting C compile features
# -- Detecting C compile features - done
# -- Detecting CXX compiler ABI info
# -- Detecting CXX compiler ABI info - done
# -- Check for working CXX compiler: /usr/bin/c++ - skipped
# -- Detecting CXX compile features
# -- Detecting CXX compile features - done
# -- Starting TileDB-SOMA superbuild.
# -- Could NOT find TileDB (missing: TileDB_DIR)
# -- Adding TileDB as an external project
# -- Could NOT find spdlog (missing: spdlog_DIR)
# -- Adding spdlog as an external project
# -- Not found clang-tidy
# -- Not found clang-format
# -- Configuring done (0.2s)
# -- Generating done (0.0s)
# -- Build files have been written to: /TileDB-SOMA/build
# [  4%] Creating directories for 'ep_tiledb'
# [  8%] Performing download step (download, verify and extract) for 'ep_tiledb'
# -- ep_tiledb download command succeeded.  See also /TileDB-SOMA/build/externals/src/ep_tiledb-stamp/ep_tiledb-download-*.log
# [ 12%] No update step for 'ep_tiledb'
# [ 16%] No patch step for 'ep_tiledb'
# [ 20%] Performing configure step for 'ep_tiledb'
# CMake Error at /TileDB-SOMA/build/externals/src/ep_tiledb-stamp/ep_tiledb-configure-Release.cmake:49 (message):
#   Command failed: 1
#
#    '/usr/bin/cmake' '-DCMAKE_INSTALL_PREFIX=/TileDB-SOMA/build/externals/install' '-DCMAKE_PREFIX_PATH=/TileDB-SOMA/build/externals/install' '-DTILEDB_S3=ON' '-DTILEDB_AZURE=ON' '-DTILEDB_GCS=OFF' '-DTILEDB_HDFS=OFF' '-DTILEDB_SERIALIZATION=ON' '-DTILEDB_WERROR=OFF' '-DTILEDB_REMOVE_DEPRECATIONS=OFF' '-DTILEDB_VERBOSE=OFF' '-DTILEDB_TESTS=OFF' '-DCMAKE_BUILD_TYPE=Release' '-DCMAKE_OSX_ARCHITECTURES=' '-DCMAKE_C_FLAGS=' '-DCMAKE_CXX_FLAGS=' '-DCMAKE_CXX_COMPILER=/usr/bin/c++' '-DCMAKE_C_COMPILER=/usr/bin/cc' '' '-DCMAKE_POSITION_INDEPENDENT_CODE=ON' '-GUnix Makefiles' '-S' '/TileDB-SOMA/build/externals/src/ep_tiledb' '-B' '/TileDB-SOMA/build/externals/src/ep_tiledb-build'
#
#   See also
#
#     /TileDB-SOMA/build/externals/src/ep_tiledb-stamp/ep_tiledb-configure-*.log
#
#
# gmake[3]: *** [CMakeFiles/ep_tiledb.dir/build.make:92: externals/src/ep_tiledb-stamp/ep_tiledb-configure] Error 1
# gmake[2]: *** [CMakeFiles/Makefile2:89: CMakeFiles/ep_tiledb.dir/all] Error 2
# gmake[1]: *** [Makefile:91: all] Error 2
# make: *** [Makefile:19: install] Error 2
```

```bash
cat /TileDB-SOMA/build/externals/src/ep_tiledb-stamp/ep_tiledb-configure-err.log
# TILEDB_DISABLE_AUTO_VCPKG is not defined. Fetch a local copy of vcpkg.
# Vcpkg commit string used: 7aeffc91033ad35cc4e2c152f213a866ec6c11ac
# Using vcpkg features: azure;serialization;s3;webp
# CMake Error at /TileDB-SOMA/build/externals/src/ep_tiledb-build/_deps/vcpkg-src/scripts/buildsystems/vcpkg.cmake:902 (message):
#   vcpkg install failed.  See logs for more information:
#   /TileDB-SOMA/build/externals/src/ep_tiledb-build/vcpkg-bootstrap.log
# Call Stack (most recent call first):
#   /usr/share/cmake-3.28/Modules/CMakeDetermineSystem.cmake:170 (include)
#   CMakeLists.txt:121 (project)
#
#
# CMake Error: CMake was unable to find a build program corresponding to "Unix Makefiles".  CMAKE_MAKE_PROGRAM is not set.  You probably need to select a different build tool.
```

```bash
tail /TileDB-SOMA/build/externals/src/ep_tiledb-build/vcpkg-bootstrap.log
# [208/216] Building CXX object CMakeFiles/vcpkg-test.dir/src/vcpkg-test/unicode.cpp.o
# [209/216] Building CXX object CMakeFiles/vcpkg-test.dir/src/vcpkg-test/xunitwriter.cpp.o
# [210/216] Building CXX object CMakeFiles/vcpkg-test.dir/src/vcpkg-test/update.cpp.o
# [211/216] Building CXX object CMakeFiles/vcpkg-test.dir/src/vcpkg-test/versionplan.cpp.o
# [212/216] Building CXX object CMakeFiles/vcpkg-test.dir/src/vcpkg-test/util.cpp.o
# [213/216] Building CXX object CMakeFiles/vcpkg-test.dir/src/vcpkg-test/plan.cpp.o
# [214/216] Building CXX object CMakeFiles/vcpkg-test.dir/src/vcpkg-test/dependencies.cpp.o
# [215/216] Building CXX object CMakeFiles/vcpkg-test.dir/src/vcpkg-test/catch.cpp.o
# [216/216] Linking CXX executable vcpkg-test
# Environment variable VCPKG_FORCE_SYSTEM_BINARIES must be set on arm, s390x, ppc64le and riscv platforms.
```

This sequence of log files is commonly useful:
- `ep_tiledb-configure-err.log`
- `vcpkg-bootstrap.log`


[vcpkg]: https://vcpkg.io/
