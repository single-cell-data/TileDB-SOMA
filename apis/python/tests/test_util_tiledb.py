import pyarrow as pa
import pytest

import tiledbsoma as soma
import tiledb


def test_stats(tmp_path, capsys: pytest.CaptureFixture[str]):
    """Make sure these exist, don't throw, and write correctly."""
    tiledb.stats_enable()
    tiledb.stats_reset()

    schema = pa.schema([("soma_joinid", pa.int64())])
    with soma.DataFrame.create(
        tmp_path.as_posix(),
        schema=schema,
        index_column_names=["soma_joinid"],
    ) as sidf:
        data = {
            "soma_joinid": [0],
        }
        sidf.write(pa.Table.from_pydict(data))

    with soma.DataFrame.open(tmp_path.as_posix()) as sidf:
        sidf.read().concat()

    tiledb.stats_dump()
    tiledb.stats_disable()
    stdout, stderr = capsys.readouterr()
    assert stdout != ""
    assert stderr == ""
    print(f"tiledbsoma_stats_dump() = {stdout}")
    soma.show_package_versions()
