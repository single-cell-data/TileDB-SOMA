import anndata as ad
import numpy as np
import pyarrow as pa
import pytest
from scipy import sparse as sp

import tiledbsoma as soma
import tiledbsoma.io as somaio


@pytest.fixture
def src_matrix(request):
    format, shape, density = request.param
    if format == "dense":
        return np.random.default_rng().standard_normal(np.prod(shape), dtype=np.float32).reshape(shape)

    return sp.random(10, 89, density=density, format=format, dtype=np.float32)


@pytest.mark.parametrize(
    "num_rows",
    [0, 1, 2, 3, 4, 10, 100, 1_000],
)
@pytest.mark.parametrize(
    "cap_nbytes",
    [1, 100, 1_000],
)
@pytest.mark.parametrize(
    "schema",
    [
        pa.schema(
            [
                ("bool", pa.bool_()),
                ("int32", pa.int32()),
                ("float64", pa.float64()),
                ("string", pa.string()),
            ],
        ),
        pa.schema(
            [
                ("bool_enum", pa.dictionary(pa.int8(), pa.bool_())),
                ("int32_enum", pa.dictionary(pa.int8(), pa.int32())),
                ("float64_enum", pa.dictionary(pa.int8(), pa.float64())),
                ("string_enum", pa.dictionary(pa.int8(), pa.string())),
            ],
        ),
    ],
)
def test_write_arrow_table(tmp_path, num_rows, cap_nbytes, schema):
    """
    Additional focus-testing for tiledbsoma.io._write_arrow_table
    """

    pydict = {}
    pydict["soma_joinid"] = list(range(num_rows))
    pydict["bool"] = [bool(e % 2) for e in range(num_rows)]
    pydict["int32"] = [(e + 1) * 10 for e in range(num_rows)]
    pydict["float64"] = [(e + 1) / 25 for e in range(num_rows)]
    pydict["string"] = [str((e + 1) / 25) for e in range(num_rows)]
    pydict["bool_enum"] = [bool(e % 2) for e in range(num_rows)]
    pydict["int32_enum"] = [((e + 1) % 24) * 10 for e in range(num_rows)]
    pydict["float64_enum"] = [((e + 1) % 24) / 25 for e in range(num_rows)]
    pydict["string_enum"] = [str(((e + 1) % 24) / 25) for e in range(num_rows)]

    tcopt = soma.TileDBCreateOptions(remote_cap_nbytes=cap_nbytes)
    twopt = soma.TileDBWriteOptions()
    uri = tmp_path.as_posix()

    table = pa.Table.from_pydict(pydict, schema=schema.insert(0, pa.field("soma_joinid", pa.int64())))
    expect_error = cap_nbytes < table[:1].nbytes  # Not enough room for even one row
    domain = [[0, max(1, len(table) - 1)]]

    with soma.DataFrame.create(uri, schema=schema, domain=domain) as sdf:
        if expect_error:
            with pytest.raises(soma.SOMAError):
                somaio.ingest._write_arrow_table(table, sdf, tcopt, twopt)
        else:
            somaio.ingest._write_arrow_table(table, sdf, tcopt, twopt)

    if not expect_error:
        with soma.DataFrame.open(uri) as sdf:
            pdf = sdf.read().concat().to_pandas()

            for column in pdf.columns.to_list():
                assert list(pdf[column]) == pydict[column]


@pytest.mark.parametrize(
    "num_rows",
    [0, 1, 2, 3, 4, 10, 100, 1_000],
)
@pytest.mark.parametrize(
    "cap_nbytes",
    [1, 100, 1_000],
)
def test_write_arrow_table_enum_to_values(tmp_path, num_rows, cap_nbytes):
    """
    Additional focus-testing for tiledbsoma.io._write_arrow_table
    """

    array_schema = pa.schema(
        [
            ("bool", pa.bool_()),
            ("int32", pa.int32()),
            ("float64", pa.float64()),
            ("string", pa.string()),
        ],
    )

    data_schema = pa.schema(
        [
            ("bool", pa.dictionary(pa.int8(), pa.bool_())),
            ("int32", pa.dictionary(pa.int8(), pa.int32())),
            ("float64", pa.dictionary(pa.int8(), pa.float64())),
            ("string", pa.dictionary(pa.int8(), pa.string())),
        ],
    )

    pydict = {}
    pydict["soma_joinid"] = list(range(num_rows))
    pydict["bool"] = [bool(e % 2) for e in range(num_rows)]
    pydict["int32"] = [((e + 1) % 24) * 10 for e in range(num_rows)]
    pydict["float64"] = [((e + 1) % 24) / 25 for e in range(num_rows)]
    pydict["string"] = [str(((e + 1) % 24) / 25) for e in range(num_rows)]

    tcopt = soma.TileDBCreateOptions(remote_cap_nbytes=cap_nbytes)
    twopt = soma.TileDBWriteOptions()
    uri = tmp_path.as_posix()

    table = pa.Table.from_pydict(pydict, schema=data_schema.insert(0, pa.field("soma_joinid", pa.int64())))
    expect_error = cap_nbytes < table[:1].nbytes  # Not enough room for even one row
    domain = [[0, max(1, len(table) - 1)]]

    with soma.DataFrame.create(uri, schema=array_schema, domain=domain) as sdf:
        if expect_error:
            with pytest.raises(soma.SOMAError):
                somaio.ingest._write_arrow_table(table, sdf, tcopt, twopt)
        else:
            somaio.ingest._write_arrow_table(table, sdf, tcopt, twopt)

    if not expect_error:
        with soma.DataFrame.open(uri) as sdf:
            pdf = sdf.read().concat().to_pandas()

            for column in pdf.columns.to_list():
                assert list(pdf[column]) == pydict[column]


def test_add_matrices(tmp_path, conftest_pbmc_small_h5ad_path):
    """Test multiple add_matrix_to_collection calls can be issued on the same soma object.

    See https://github.com/single-cell-data/TileDB-SOMA/issues/1565."""
    # Create a soma object from an anndata object
    soma_uri = soma.io.from_h5ad(
        tmp_path.as_posix(),
        input_path=conftest_pbmc_small_h5ad_path,
        measurement_name="RNA",
    )

    # Synthesize some new data to be written into two matrices within the soma object (ensuring it's different from the
    # original data, so that writes must be performed)
    h5ad = ad.read_h5ad(conftest_pbmc_small_h5ad_path)
    new_X_pca = h5ad.obsm["X_pca"] * 2
    new_PCs = h5ad.varm["PCs"] * 2

    with soma.open(soma_uri, "w") as exp_w:
        # Write multiple new matrices into the soma object, attempting to repro https://github.com/single-cell-data/TileDB-SOMA/issues/1565
        soma.io.add_matrix_to_collection(
            exp=exp_w,
            measurement_name="RNA",
            collection_name="obsm",
            matrix_name="logcounts_pca",
            matrix_data=new_X_pca,
        )

        soma.io.add_matrix_to_collection(
            exp=exp_w,
            measurement_name="RNA",
            collection_name="varm",
            matrix_name="logcounts_pcs",
            matrix_data=new_PCs,
        )
