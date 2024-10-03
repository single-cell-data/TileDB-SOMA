# This file enables building the libtiledbsoma external package as part of the
# tiledbsoma module install process.
#
# Local non-editable install:
#   `pip install .`
#
# Local editable install:
#   `pip install -e .`
#
# Install from PyPI
#   `pip install tiledbsoma`
#
# Based on ideas from https://github.com/pybind/cmake_example
# The `bld` script here is reused for pip install, CI, and local builds.

import ctypes
import os
import pathlib
import shutil
import subprocess
import sys
from typing import Optional

import setuptools.command.bdist_wheel
import setuptools.command.build_ext

try:
    from pybind11.setup_helpers import Pybind11Extension

except ImportError:
    # Explanation:
    # https://pybind11.readthedocs.io/en/stable/compiling.html#classic-setup-requires
    # This works around a catch-22 where pybind11 is a requirement to load setup.py, yet pip cannot
    # read&install our requirements without loading setup.py.
    from setuptools import Extension as Pybind11Extension

this_dir = pathlib.Path(__file__).parent.absolute()
sys.path.insert(0, str(this_dir))

import version  # noqa E402

# tiledb_dir and tiledbsoma_dir may be specified by either environment variable
# or command-line argument. If both are provided, the latter wins.

tiledb_dir: Optional[pathlib.Path] = None
tiledbsoma_dir: Optional[pathlib.Path] = None
no_tiledb_dep: bool = False

args = sys.argv[:]
for arg in args:
    start, eq, last = arg.partition("=")
    if (start, eq) == ("--tiledb", "="):
        tiledb_dir = pathlib.Path(last)
        sys.argv.remove(arg)
    if (start, eq) == ("--libtiledbsoma", "="):
        tiledbsoma_dir = pathlib.Path(last)
        sys.argv.remove(arg)
    if (start, eq) == ("--no-tiledb-deprecated", "="):
        no_tiledb_dep = True
        sys.argv.remove(arg)


tiledb_dir = os.environ.get("TILEDB_PATH", tiledb_dir)
tiledb_given = tiledb_dir is not None
tiledbsoma_dir = os.environ.get("TILEDBSOMA_PATH", tiledbsoma_dir)

if tiledbsoma_dir is not None and tiledb_dir is None:
    raise ValueError(
        "If TILEDBSOMA_PATH is set, then TILEDB_PATH must "
        "also be set. TILEDB_PATH must be set to the location of "
        "the TileDB shared object library linked to the "
        "TileDB-SOMA shared object library"
    )

if tiledbsoma_dir is None:
    scripts_dir = this_dir / "dist_links" / "scripts"
    scripts_dir = scripts_dir.resolve()

    tiledbsoma_dir = scripts_dir.parent / "dist"
else:
    tiledbsoma_dir = pathlib.Path(tiledbsoma_dir)

if tiledb_dir is None:
    # tiledb_dir = pathlib.Path(this_dir)
    scripts_dir = this_dir / "dist_links" / "scripts"
    scripts_dir = scripts_dir.resolve()

    tiledb_dir = scripts_dir.parent / "dist"
else:
    tiledb_dir = pathlib.Path(tiledb_dir)


def get_libtiledbsoma_library_name():
    """
    :return: List of TileDB shared library names.
    """
    if os.name == "posix":
        if sys.platform == "darwin":
            return "libtiledbsoma.dylib"
        else:
            return "libtiledbsoma.so"
    elif os.name == "nt":
        return "tiledbsoma.dll"
    else:
        raise RuntimeError(f"Unsupported OS name {os.name}")


def find_libtiledbsoma_full_path_on_linux(lib_name):
    # https://stackoverflow.com/questions/35682600/get-absolute-path-of-shared-library-in-python
    class LINKMAP(ctypes.Structure):
        _fields_ = [("l_addr", ctypes.c_void_p), ("l_name", ctypes.c_char_p)]

    libdl = ctypes.CDLL(lib_name)
    dlinfo = libdl.dlinfo
    dlinfo.argtypes = ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p
    dlinfo.restype = ctypes.c_int

    libc = ctypes.CDLL(lib_name)
    lmptr = ctypes.c_void_p()
    dlinfo(libc._handle, 2, ctypes.byref(lmptr))

    return ctypes.cast(lmptr, ctypes.POINTER(LINKMAP)).contents.l_name.decode()


def libtiledbsoma_exists():
    """
    Returns the path to the TileDB-SOMA library, if it exists.
    :return: The path to the TileDB-SOMA library, or None.
    """
    # Check if TileDB-SOMA is installed in user given path
    dist_dirs = [tiledbsoma_dir / "lib"]
    if sys.platform.startswith("linux"):
        dist_dirs.append(tiledbsoma_dir / "lib64")
        dist_dirs.append(tiledbsoma_dir / "lib" / "x86_64-linux-gnu")
    elif os.name == "nt":
        dist_dirs.append(tiledbsoma_dir / "bin")

    for lib_dir in dist_dirs:
        full_lib_path = lib_dir / get_libtiledbsoma_library_name()
        print(f"Checking: {full_lib_path} exists: {full_lib_path.exists()}")
        if full_lib_path.exists():
            return lib_dir

    # Check to see if TileDB-SOMA is globally installed.
    lib_name = get_libtiledbsoma_library_name()

    try:
        # Note: This is a relative path on Linux
        # https://bugs.python.org/issue21042
        if os.name == "posix" and sys.platform != "darwin":
            path = find_libtiledbsoma_full_path_on_linux(lib_name)
        else:
            path = ctypes.CDLL(lib_name)
        print(f"Found globally installed {path}")
        return pathlib.Path(path).parents[0]
    except Exception as e:
        print(e)
        return None


def find_or_build_package_data(setuptools_cmd):
    # check if libtiledbsoma is installed
    lib_dir = libtiledbsoma_exists()

    # if not then build from source
    if lib_dir is None:
        # Note: The GitHub build process uses the contents of `bld` as a key
        # to cache the native binaries. Using non-default options here will
        # cause that cache to fall out of sync.
        #
        # See `.github/workflows/python-ci-single.yml` for configuration.
        if os.name == "nt":
            bld_command = ["pwsh.exe", "./bld.ps1"]
            if tiledb_dir is not None:
                bld_command.append(f"TileDBLocation={tiledb_dir}")
            if no_tiledb_dep:
                bld_command.append("RemoveTileDBDeprecated=ON")

        else:
            bld_command = ["./bld"]
            if tiledb_dir is not None:
                bld_command.append(f"--tiledb={tiledb_dir}")
            if no_tiledb_dep:
                bld_command.append("--no-tiledb-deprecated=true")

        subprocess.run(bld_command, cwd=scripts_dir, check=True)

        lib_dir = libtiledbsoma_exists()
        assert lib_dir, "error when building libtiledbsoma from source"

        # If we are building from source, then we are likely building wheels.
        # Copy both the tiledbsoma and tiledb shared objects into the
        # package dir so they can be found by package_data
        package_data = []
        src_dir = this_dir / "src" / "tiledbsoma"
        for f in lib_dir.glob("*tiledbsoma.*"):
            if f.suffix != ".a":  # skip static library
                print(f"  copying file {f} to {src_dir}")
                shutil.copy(f, src_dir)
                package_data.append(f.name)
        assert package_data, f"tiledbsoma artifacts absent from {lib_dir}"

        if not tiledb_given:
            for f in lib_dir.glob("*tiledb.*"):
                if f.suffix != ".a":  # skip static library
                    print(f"  copying file {f} to {src_dir}")
                    shutil.copy(f, src_dir)
                    package_data.append(f.name)
                assert package_data, f"tiledb artifacts absent from {lib_dir}"

        # Install shared libraries inside the Python module via package_data.
        print(f"  adding to package_data: {package_data}")
        setuptools_cmd.distribution.package_data["tiledbsoma"] = package_data


class build_ext(setuptools.command.build_ext.build_ext):
    def run(self):
        find_or_build_package_data(self)
        super().run()


class bdist_wheel(setuptools.command.bdist_wheel.bdist_wheel):
    def run(self):
        find_or_build_package_data(self)
        super().run()


INC_DIRS = [
    "dist_links/libtiledbsoma/include",
    "dist_links/libtiledbsoma/external/include",
    "dist_links/libtiledbsoma/src/external/include",
    "../../build/externals/install/include",
    str(tiledbsoma_dir / "include"),
    str(tiledbsoma_dir / "include/tiledbsoma"),
    str(tiledbsoma_dir.parent / "build/externals/install/include"),
    str(tiledbsoma_dir / "include"),
    str(tiledb_dir / "include"),
]

LIB_DIRS = [
    str(tiledbsoma_dir / "lib"),
    str(tiledb_dir / "lib"),
]

CXX_FLAGS = []

if os.name != "nt":
    CXX_FLAGS.append(f'-Wl,-rpath,{str(tiledbsoma_dir / "lib")}')
    CXX_FLAGS.append(f'-Wl,-rpath,{str(tiledb_dir / "lib")}')

if sys.platform == "darwin":
    CXX_FLAGS.append("-mmacosx-version-min=11.0")

if os.name == "posix" and sys.platform != "darwin":
    LIB_DIRS.append(str(tiledbsoma_dir / "lib" / "x86_64-linux-gnu"))
    LIB_DIRS.append(str(tiledbsoma_dir / "lib64"))
    LIB_DIRS.append(str(tiledb_dir / "lib" / "x86_64-linux-gnu"))
    LIB_DIRS.append(str(tiledb_dir / "lib64"))
    CXX_FLAGS.append(f'-Wl,-rpath,{str(tiledbsoma_dir / "lib" / "x86_64-linux-gnu")}')
    CXX_FLAGS.append(f'-Wl,-rpath,{str(tiledbsoma_dir / "lib64")}')
    CXX_FLAGS.append(f'-Wl,-rpath,{str(tiledb_dir / "lib" / "x86_64-linux-gnu")}')
    CXX_FLAGS.append(f'-Wl,-rpath,{str(tiledb_dir / "lib64")}')

# ----------------------------------------------------------------
# Don't use `if __name__ == "__main__":` as the `python_requires` must
# be at top level, outside any if-block
# https://github.com/pypa/cibuildwheel/blob/7c4bbf8cb31d856a0fe547faf8edf165cd48ce74/cibuildwheel/projectfiles.py#L41-L46

setuptools.setup(
    name="tiledbsoma",
    description="Python API for efficient storage and retrieval of single-cell data using TileDB",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="TileDB, Inc.",
    author_email="help@tiledb.io",
    maintainer="TileDB, Inc.",
    maintainer_email="help@tiledb.io",
    url="https://github.com/single-cell-data/TileDB-SOMA/tree/main/apis/python",
    license="MIT",
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Information Technology",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Operating System :: Unix",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages("src"),
    # This next is necessary to avoid cibuildwheel thinking we want a python-only wheel:
    ext_modules=[
        Pybind11Extension(
            "tiledbsoma.pytiledbsoma",
            [
                "src/tiledbsoma/common.cc",
                "src/tiledbsoma/reindexer.cc",
                "src/tiledbsoma/query_condition.cc",
                "src/tiledbsoma/soma_context.cc",
                "src/tiledbsoma/soma_array.cc",
                "src/tiledbsoma/soma_object.cc",
                "src/tiledbsoma/soma_dataframe.cc",
                "src/tiledbsoma/soma_point_cloud_dataframe.cc",
                "src/tiledbsoma/soma_dense_ndarray.cc",
                "src/tiledbsoma/soma_sparse_ndarray.cc",
                "src/tiledbsoma/soma_group.cc",
                "src/tiledbsoma/soma_collection.cc",
                "src/tiledbsoma/pytiledbsoma.cc",
            ],
            include_dirs=INC_DIRS,
            library_dirs=LIB_DIRS,
            libraries=["tiledbsoma"] + (["tiledb"] if os.name == "nt" else []),
            extra_link_args=CXX_FLAGS,
            extra_compile_args=["-std=c++17" if os.name != "nt" else "/std:c++17"]
            + CXX_FLAGS,
            language="c++",
        )
    ],
    zip_safe=False,
    setup_requires=["pybind11"],
    install_requires=[
        # Tracked in https://github.com/single-cell-data/TileDB-SOMA/issues/1785
        "anndata != 0.10.0",
        "attrs>=22.2",
        "numba>=0.58.0",
        "numpy<2.0",
        "pandas",
        "pyarrow",
        "scanpy>=1.9.2",
        "scipy",
        # Note: the somacore version is in .pre-commit-config.yaml too
        "somacore==1.0.20",
        "tiledb~=0.32.0",
        "typing-extensions",  # Note "-" even though `import typing_extensions`
    ],
    extras_require={
        "dev": open("requirements_dev.txt").read(),
    },
    python_requires=">=3.9",
    cmdclass={"build_ext": build_ext, "bdist_wheel": bdist_wheel},
    version=version.get_version(),
)
