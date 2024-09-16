import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import tiledbsoma as soma
import tiledb
import tempfile

with tempfile.TemporaryDirectory() as tmpdirname:
    print('created temporary directory', tmpdirname)

    b = pq.read_table('/mnt/c/Users/xanthos/Downloads/Xenium_V1_FFPE_TgCRND8_2_5_months_outs/cells.parquet')

    print(soma.GeometryDataFrame.create(tmpdirname, schema=b.schema))

    a = pq.read_table('/mnt/c/Users/xanthos/Downloads/Xenium_V1_FFPE_TgCRND8_2_5_months_outs/cell_boundaries.parquet')
    b = pq.read_table('/mnt/c/Users/xanthos/Downloads/Xenium_V1_FFPE_TgCRND8_2_5_months_outs/cells.parquet')
    obj = soma.GeometryDataFrame.open(tmpdirname)

    print(soma.open(tmpdirname))