

cmake_bin <- Sys.which("cmake")
if (cmake_bin == "" && file.exists("/Applications/CMake.app/Contents/bin/cmake"))
    cmake_bin <- "/Applications/CMake.app/Contents/bin/cmake"

if (cmake_bin == "") {
    stop("No 'cmake' binary found", call. = FALSE)
}

git_bin <- Sys.which("git")
if (git_bin == "") {
    stop("No 'git' binary found", call. = FALSE)
}

lines <- readLines("tools/build_libtiledbsoma.sh.in")
lines <- gsub("@cmake@", cmake_bin, lines)
writeLines(lines, "tools/build_libtiledbsoma.sh")
Sys.chmod("tools/build_libtiledbsoma.sh", mode="0755")
