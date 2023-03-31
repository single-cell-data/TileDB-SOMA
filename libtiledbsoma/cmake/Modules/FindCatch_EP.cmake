# catch2 has a set-but-unused-variable warning as of MacOS 13.3.
# So we need to not use -Werror, as that's third-party code out of our control.
#
# See also:
# https://discourse.cmake.org/t/how-to-turn-off-warning-flags-for-project-added-by-fetchcontent-declare/2461

Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.0.1
)

get_property(
    compile_options
    DIRECTORY
    PROPERTY COMPILE_OPTIONS
)

set_property(
    DIRECTORY
    APPEND
    PROPERTY COMPILE_OPTIONS -Wno-error=unused-but-set-variable
)

FetchContent_MakeAvailable(Catch2)

set_property(
    DIRECTORY
    PROPERTY COMPILE_OPTIONS ${compile_options}
)

unset(compile_options)
