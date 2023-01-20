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
        "numpy",
        "pandas",
        "pyarrow >= 9.0.0",
        "scanpy",
        "scipy",
        "somacore==0.0.0a5",
        "tiledb>=0.20.0",
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
    cmdclass={"build_ext": build_ext, "bdist_wheel": bdist_wheel},
    version=version.getVersion(),
)
