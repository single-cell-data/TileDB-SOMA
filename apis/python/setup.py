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

MODULE_NAME = "tiledbsoma"
EXT_NAME = "tiledbsoma.libtiledbsoma"


def find_or_build(setuptools_cmd):
    # TODO: support windows
    if sys.platform.startswith("win32"):
        return

    # Setup paths
    python_dir = os.path.abspath(os.path.dirname(__file__))
    src_dir = f"{python_dir}/src/{MODULE_NAME}"
    scripts_dir = f"{python_dir}/../../scripts"
    lib_dir = f"{python_dir}/../../dist/lib"

    # Call the build script if the install library directory does not exist
    if not os.path.exists(lib_dir):
        subprocess.check_call([f"{scripts_dir}/bld"])

    # Copy native libs into the package dir so they can be found by package_data
    package_data = []
    for obj in [os.path.join(lib_dir, f) for f in os.listdir(lib_dir)]:
        print(f"  copying file {obj} to {src_dir}")
        shutil.copy(obj, src_dir)
        package_data.append(os.path.basename(obj))

    # Install shared libraries inside the Python module via package_data.
    print(f"  adding to package_data: {package_data}")
    setuptools_cmd.distribution.package_data.update({MODULE_NAME: package_data})


def get_ext_modules():
    # Workaround windows issue compiling an extension with no source files
    if sys.platform.startswith("win32"):
        return []
    else:
        return [CMakeExtension(EXT_NAME)]


class CMakeExtension(Extension):
    def __init__(self, name):
        Extension.__init__(self, name, sources=[])


class BuildExtCmd(build_ext):
    def run(self):
        find_or_build(self)
        super().build_extensions()


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
        setup_requires=[
            "setuptools>=65.3.0",
            "setuptools_scm>=1.5.4",
            "wheel>=0.30",
            "setuptools_scm_git_archive",
        ],
        install_requires=[
            "anndata",
            "pandas",
            "pyarrow",
            "scanpy",
            "scipy",
            "tiledb>=0.17.0",
            "pybind11>=2.10.0",
            "pytest>=7.1.3",
        ],
        python_requires=">=3.7",
        ext_modules=get_ext_modules(),
        cmdclass={
            "build_ext": BuildExtCmd,
            "bdist_egg": BdistEggCmd,
            "bdist_wheel": BdistWheelCmd,
        },
    )
