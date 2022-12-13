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

import os
import shutil
import subprocess
import sys

from setuptools import Extension, find_packages, setup
from setuptools.command.bdist_egg import bdist_egg
from setuptools.command.build_ext import build_ext
from wheel.bdist_wheel import bdist_wheel

sys.path.insert(0, os.path.dirname(__file__))
import version  # noqa E402

MODULE_NAME = "tiledbsoma"
EXT_NAME = "tiledbsoma.libtiledbsoma"


def find_or_build(setuptools_cmd):
    # Setup paths
    python_dir = os.path.abspath(os.path.dirname(__file__))
    src_dir = f"{python_dir}/src/{MODULE_NAME}"
    if os.path.islink(os.path.join(python_dir, "dist_links/scripts")):
        # in git source tree
        scripts_dir = f"{python_dir}/../../scripts"
        lib_dir = f"{python_dir}/../../dist/lib"
    else:
        # in extracted sdist, with libtiledbsoma copied into dist_links/
        scripts_dir = f"{python_dir}/dist_links/scripts"
        lib_dir = f"{python_dir}/dist_links/dist/lib"

    # Call the build script if the install library directory does not exist
    if not os.path.exists(lib_dir):
        subprocess.run("bash bld", cwd=scripts_dir, shell=True)

    # Copy native libs into the package dir so they can be found by package_data
    package_data = []
    for obj in [os.path.join(lib_dir, f) for f in os.listdir(lib_dir)]:
        # skip static library
        if not obj.endswith(".a"):
            print(f"  copying file {obj} to {src_dir}")
            shutil.copy(obj, src_dir)
            package_data.append(os.path.basename(obj))

    # Install shared libraries inside the Python module via package_data.
    print(f"  adding to package_data: {package_data}")
    setuptools_cmd.distribution.package_data.update({MODULE_NAME: package_data})


def get_ext_modules():
    return [CMakeExtension(EXT_NAME)]


class CMakeExtension(Extension):
    def __init__(self, name):
        Extension.__init__(self, name, sources=[])


class BuildExtCmd(build_ext):
    def run(self):
        find_or_build(self)


class BdistEggCmd(bdist_egg):
    def run(self):
        find_or_build(self)
        bdist_egg.run(self)


class BdistWheelCmd(bdist_wheel):
    def run(self):
        find_or_build(self)
        bdist_wheel.run(self)


if __name__ == "__main__":
    setup(
        name=MODULE_NAME,
        description="Python API for efficient storage and retrieval of single-cell data using TileDB",
        author="TileDB, Inc.",
        author_email="help@tiledb.io",
        maintainer="TileDB, Inc.",
        maintainer_email="help@tiledb.io",
        url="https://github.com/single-cell-data/TileDB-SOMA/apis/python",
        license="MIT",
        classifiers=[
            "Intended Audience :: Developers",
            "Intended Audience :: Information Technology",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Operating System :: Unix",
            "Operating System :: POSIX :: Linux",
            "Operating System :: MacOS :: MacOS X",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
        ],
        package_dir={"": "src"},
        packages=find_packages("src"),
        zip_safe=False,
        install_requires=[
            "anndata",
            "pandas",
            "pyarrow",
            "scanpy",
            "scipy",
            "tiledb>=0.19.0",
        ],
        python_requires=">=3.7",
        ext_modules=get_ext_modules(),
        cmdclass={
            "build_ext": BuildExtCmd,
            "bdist_egg": BdistEggCmd,
            "bdist_wheel": BdistWheelCmd,
        },
        version=version.getVersion(),
    )
