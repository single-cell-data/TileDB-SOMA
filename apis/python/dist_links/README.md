# dist_links

This directory contains symlinks to sibling directories in this repo that are used by our Python wrappers. They trigger those directories' inclusion in our [sdist tarballs](https://docs.python.org/3/distutils/sourcedist.html), in which all the necessary sources should be nested inside the Python package directory. The links make it appear this is the case to [`setup.py`](https://github.com/single-cell-data/TileDB-SOMA/blob/main/apis/python/setup.py) without actually reorganizing our repo.
