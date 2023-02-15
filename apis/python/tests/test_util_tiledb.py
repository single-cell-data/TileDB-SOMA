import pyarrow as pa
import pytest

import tiledbsoma as soma


def test_stats(tmp_path, capfd: pytest.CaptureFixture[str]):
    """Make sure these exist, don't throw, and write correctly."""
    soma.tiledbsoma_stats_enable()
    soma.tiledbsoma_stats_reset()

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

    soma.tiledbsoma_stats_dump()
    soma.tiledbsoma_stats_disable()
    stdout, stderr = capfd.readouterr()
    assert stdout != ""
    assert stderr == ""
