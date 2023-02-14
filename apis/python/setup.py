# This file is enables building the libtiledbsoma external package as part of the
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

# type: ignore
import ctypes
import os
import pathlib
import shutil
import subprocess
import sys

import setuptools
import setuptools.command.build_ext
import wheel.bdist_wheel

this_dir = pathlib.Path(__file__).parent.absolute()
sys.path.insert(0, str(this_dir))
import version  # noqa E402


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
    Returns the path to the globally installed TileDB-SOMA library, if it exists.
    :return: The path to the TileDB-SOMA library, or None.
    """
    # Check to see if TileDB-SOMA is globally installed.
    lib_name = get_libtiledbsoma_library_name()

    try:
        # Note: This is a relative path on Linux
        # https://bugs.python.org/issue21042
        if os.name == "posix" and sys.platform != "darwin":
            path = find_libtiledbsoma_full_path_on_linux(lib_name)
        else:
            path = ctypes.CDLL(lib_name)
        return pathlib.Path(path).parents[0]
    except Exception as e:
        print(e)
        return None


def find_or_build_package_data(setuptools_cmd):
    global libtiledbsoma_dir

    # Set up paths
    scripts_dir = this_dir / "dist_links" / "scripts"

    if scripts_dir.is_symlink():
        # in git source tree
        libtiledbsoma_dir = this_dir.parent.parent
    else:
        # in extracted sdist, with libtiledbsoma copied into dist_links/
        libtiledbsoma_dir = this_dir / "dist_links"

    # Call the build script if the install library directory does not exist
    lib_dir = libtiledbsoma_dir / "dist" / "lib"

    if not lib_dir.exists():
        # check if libtilesoma is globally installed
        global_libtiledbsoma_path = libtiledbsoma_exists()
        if global_libtiledbsoma_path is not None:
            lib_dir = global_libtiledbsoma_path
        else:
            # If not then build from source

            # Note: The GitHub build process uses the contents of `bld` as a key
            # to cache the native binaries. Using non-default options here will
            # cause that cache to fall out of sync.
            #
            # See `.github/workflows/python-ci-single.yml` for configuration.
            subprocess.run(["./bld"], cwd=scripts_dir)

    # Copy native libs into the package dir so they can be found by package_data
    package_data = []
    src_dir = this_dir / "src" / "tiledbsoma"
    for f in lib_dir.glob("*tiledbsoma*"):
        if f.suffix != ".a":  # skip static library
            print(f"  copying file {f} to {src_dir}")
            shutil.copy(f, src_dir)
            package_data.append(f.name)

    # Install shared libraries inside the Python module via package_data.
    print(f"  adding to package_data: {package_data}")
    setuptools_cmd.distribution.package_data["tiledbsoma"] = package_data


class build_ext(setuptools.command.build_ext.build_ext):
    def run(self):
        find_or_build_package_data(self)


class bdist_wheel(wheel.bdist_wheel.bdist_wheel):
    def run(self):
        find_or_build_package_data(self)
        super().run()


# ----------------------------------------------------------------
# Don't use `if __name__ == "__main__":` as the `python_requires` must
# be at top level, outside any if-block
# https://github.com/pypa/cibuildwheel/blob/7c4bbf8cb31d856a0fe547faf8edf165cd48ce74/cibuildwheel/projectfiles.py#L41-L46
setuptools.setup(
    name="tiledbsoma",
    description="Python API for efficient storage and retrieval of single-cell data using TileDB",
    long_description=open("README.md").read(),
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
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages("src"),
    # This next is necessary to avoid cibuildwheel thinking we want a python-only wheel:
    ext_modules=[setuptools.Extension("tiledbsoma.libtiledbsoma", sources=[])],
    zip_safe=False,
    install_requires=[
        "anndata",
        "attrs>=22.1",
        "matplotlib<3.7",  # https://github.com/scverse/scanpy/issues/2411
        "numpy",
        "pandas",
        "pyarrow >= 9.0.0",
        "scanpy",
        "scipy",
        "somacore==0.0.0a13",
        "tiledb==0.20.*",
        "typing-extensions",  # Note "-" even though `import typing_extensions`
    ],
    extras_require={
        "dev": [
            "black",
            "ruff",
            "pytest",
            "typeguard",
        ]
    },
    python_requires=">=3.7",
    cmdclass={"build_ext": build_ext, "bdist_wheel": bdist_wheel},
    version=version.getVersion(),
)
