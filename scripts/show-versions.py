#!/usr/bin/env python

import sys

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import scanpy as sc
import scipy as sp
import tiledb

import tiledbsoma

print("tiledbsoma.__version__   ", tiledbsoma.__version__)
print("tiledb.__version__       ", tiledb.__version__)
print("core version             ", ".".join(map(str, tiledb.libtiledb.version())))
print("anndata.__version__  (ad)", ad.__version__)
print("numpy.__version__    (np)", np.__version__)
print("pandas.__version__   (pd)", pd.__version__)
print("pyarrow.__version__  (pa)", pa.__version__)
print("scanpy.__version__   (sc)", sc.__version__)
print("scipy.__version__    (sp)", sp.__version__)
v = sys.version_info
print("python__version__        ", ".".join(map(str, (v.major, v.minor, v.micro))))
