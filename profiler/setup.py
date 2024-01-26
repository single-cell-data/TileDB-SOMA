from setuptools import find_packages, setup

setup(
    name="profiler",
    version="1.0",
    packages=find_packages("src"),
    package_dir={"": "src"},
    requires=["gitpython", "psutil", "tiledbsoma", "cellxgene_census"],
)
