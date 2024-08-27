

cmake_bin <- Sys.which("cmake")
if (cmake_bin == "" && file.exists("/Applications/CMake.app/Contents/bin/cmake"))
    cmake_bin <- "/Applications/CMake.app/Contents/bin/cmake"

if (cmake_bin == "") {
    message("\n*** No 'cmake' binary found.\n*** Please install 'cmake' to build from source.\n")
    q(save = "no", status = 1)
}

git_bin <- Sys.which("git")
if (git_bin == "") {
    message("\n*** No 'git' binary found.\n*** Please install 'git' to build from source.\n")
    q(save = "no", status = 1)
}

lines <- readLines("tools/build_libtiledbsoma.sh.in")
lines <- gsub("@cmake@", cmake_bin, lines)
writeLines(lines, "tools/build_libtiledbsoma.sh")
Sys.chmod("tools/build_libtiledbsoma.sh", mode="0755")
