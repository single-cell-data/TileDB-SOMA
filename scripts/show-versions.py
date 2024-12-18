#!/usr/bin/env python

import sys

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import scanpy as sc
import scipy as sp

try:
    import spatialdata as sd

    HAS_SPATIALDATA = True
except ImportError:

    HAS_SPATIALDATA = False


import tiledbsoma

print("tiledbsoma.__version__   ", tiledbsoma.__version__)
print("core version             ", tiledbsoma.get_libtiledbsoma_core_version())
print("anndata.__version__  (ad)", ad.__version__)
print("numpy.__version__    (np)", np.__version__)
print("pandas.__version__   (pd)", pd.__version__)
print("pyarrow.__version__  (pa)", pa.__version__)
print("scanpy.__version__   (sc)", sc.__version__)
print("scipy.__version__    (sp)", sp.__version__)
if HAS_SPATIALDATA:
    print("spatialdata.__version__ (sd)", sd.__version__)
else:
    print("spatialdata not installed")
v = sys.version_info
print("python__version__        ", ".".join(map(str, (v.major, v.minor, v.micro))))
