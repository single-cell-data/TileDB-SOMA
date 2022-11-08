#!/usr/bin/env python

import sys

import tiledb
import tiledbsoma
import tiledbsoma.io

t = tiledbsoma

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import scanpy as sc
import scipy as sp

print("tiledbsoma.__version__   ", tiledbsoma.__version__)
print("tiledb.__version__       ", tiledb.__version__)
print(
    "core version             ",
    ".".join(str(ijk) for ijk in list(tiledb.libtiledb.version())),
)
print("anndata.__version__  (ad)", ad.__version__)
print("numpy.__version__    (np)", np.__version__)
print("pandas.__version__   (pd)", pd.__version__)
print("pyarrow.__version__  (pa)", pa.__version__)
print("scanpy.__version__   (sc)", sc.__version__)
print("scipy.__version__    (sp)", sp.__version__)
print(
    "python__version__        ",
    ".".join(
        [
            str(e)
            for e in [
                sys.version_info.major,
                sys.version_info.minor,
                sys.version_info.micro,
            ]
        ]
    ),
)
