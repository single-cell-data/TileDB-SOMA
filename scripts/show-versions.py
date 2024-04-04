#!/usr/bin/env python

import sys

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import scanpy as sc
import scipy as sp

import tiledbsoma
import tiledb

print("tiledbsoma.__version__   ", tiledbsoma.__version__)
print("tiledb.version()         ", ".".join(str(e) for e in tiledb.version()))
print("core version             ", ".".join(map(str, tiledb.libtiledb.version())))
print("anndata.__version__  (ad)", ad.__version__)
print("numpy.__version__    (np)", np.__version__)
print("pandas.__version__   (pd)", pd.__version__)
print("pyarrow.__version__  (pa)", pa.__version__)
print("scanpy.__version__   (sc)", sc.__version__)
print("scipy.__version__    (sp)", sp.__version__)
v = sys.version_info
print("python__version__        ", ".".join(map(str, (v.major, v.minor, v.micro))))
