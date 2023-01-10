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

import pathlib
import shutil
import subprocess
import sys

import setuptools
import wheel.bdist_wheel
from setuptools.command.bdist_egg import bdist_egg
from setuptools.command.build_ext import build_ext

this_dir = pathlib.Path(__file__).parent.absolute()
sys.path.insert(0, str(this_dir))
import version  # noqa E402

MODULE_NAME = "tiledbsoma"
EXT_NAME = "tiledbsoma.libtiledbsoma"


def find_or_build_package_data(setuptools_cmd):
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
        subprocess.run("bash bld", cwd=scripts_dir, shell=True)

    # Copy native libs into the package dir so they can be found by package_data
    package_data = []
    src_dir = this_dir / "src" / "tiledbsoma"
    for f in lib_dir.glob("*"):
        if f.suffix != ".a":  # skip static library
            print(f"  copying file {f} to {src_dir}")
            shutil.copy(f, src_dir)
            package_data.append(f.name)

    # Install shared libraries inside the Python module via package_data.
    print(f"  adding to package_data: {package_data}")
    setuptools_cmd.distribution.package_data.update({MODULE_NAME: package_data})

    return package_data


def get_ext_modules():
    return [cmake_extension(EXT_NAME)]


class cmake_extension(setuptools.Extension):
    def __init__(self, name):
        setuptools.Extension.__init__(self, name, sources=[])


class build_ext_cmd(build_ext):
    def run(self):
        find_or_build_package_data(self)


class bdist_egg_cmd(bdist_egg):
    def run(self):
        find_or_build_package_data(self)
        bdist_egg.run(self)


class bdist_wheel_cmd(wheel.bdist_wheel.bdist_wheel):
    def run(self):
        package_data = find_or_build_package_data(self)
        # Install shared libraries inside the Python module via package_data
        print(f"  adding to package_data: {package_data}")
        self.distribution.package_data["tiledbsoma"] = package_data
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
    url="https://github.com/single-cell-data/TileDB-SOMA/apis/python",
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
        "numpy",
        "pandas",
        "pyarrow",
        "scanpy",
        "scipy",
        "tiledb>=0.19.0",
        "typing-extensions",  # Note "-" even though `import typing_extensions`
    ],
    extras_require={
        "dev": [
            "black",
            "flake8-bugbear",
            "isort",
            "pytest",
            "typeguard",
        ]
    },
    python_requires=">=3.7",
    cmdclass={
        "build_ext": build_ext_cmd,
        "bdist_egg": bdist_egg_cmd,
        "bdist_wheel": bdist_wheel_cmd,
    },
    version=version.getVersion(),
)
