#
# CheckAVX2Support.cmake
#
# Licensed under the MIT License.
# Copyright (c) TileDB, Inc.
#
# This file defines a function to detect toolchain support for AVX2.
#

include(CheckCXXSourceRuns)
include(CMakePushCheckState)

#
# Tries to build and run an AVX2 program with the given compiler flag.
# If successful, sets cache variable HAVE_AVX2 to 1.
#
function (CheckAVX2Flag FLAG)
  cmake_push_check_state()
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${FLAG}")
  unset(HAVE_AVX2 CACHE)
  check_cxx_source_runs("
    #include <immintrin.h>
    int main() {
      __m256i packed = _mm256_set_epi32(-1, -2, -3, -4, -5, -6, -7, -8);
      __m256i absolute_values = _mm256_abs_epi32(packed);
      return 0;
    }"
    HAVE_AVX2
  )
  cmake_pop_check_state()
endfunction()

#
# Determines if AVX2 is available.
#
# This function sets two variables in the cache:
#    COMPILER_SUPPORTS_AVX2 - Set to true if the compiler supports AVX2.
#    COMPILER_AVX2_FLAG - Set to the appropriate flag to enable AVX2.
#
function (CheckAVX2Support)
  # If defined to a false value other than "", return without checking for avx2 support
  if (DEFINED COMPILER_SUPPORTS_AVX2 AND
    NOT COMPILER_SUPPORTS_AVX2 STREQUAL "" AND
    NOT COMPILER_SUPPORTS_AVX2)
      message("AVX2 compiler support disabled by COMPILER_SUPPORTS_AVX2=${COMPILER_SUPPORTS_AVX2}")
      return()
  endif()

  if (MSVC)
    CheckAVX2Flag(/arch:AVX2)
    if (HAVE_AVX2)
      set(COMPILER_AVX2_FLAG "/arch:AVX2" CACHE STRING "Compiler flag for AVX2 support.")
    endif()
  else()
    CheckAVX2Flag(-mavx2)
    if (HAVE_AVX2)
      set(COMPILER_AVX2_FLAG "-mavx2" CACHE STRING "Compiler flag for AVX2 support.")
    endif()
  endif()

  set(COMPILER_SUPPORTS_AVX2 "${HAVE_AVX2}" CACHE BOOL "True if the compiler supports AVX2." FORCE)
  unset(HAVE_AVX2 CACHE)
endfunction()
