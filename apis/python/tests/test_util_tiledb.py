import pyarrow as pa
import pytest

import tiledbsoma


def test_stats(tmp_path, capfd: pytest.CaptureFixture[str]):
    """Make sure these exist, don't throw, and write correctly."""
    tiledbsoma.stats_enable()
    tiledbsoma.stats_reset()

    schema = pa.schema([("soma_joinid", pa.int64())])
    with tiledbsoma.DataFrame.create(
        tmp_path.as_posix(),
        schema=schema,
        index_column_names=["soma_joinid"],
    ) as sidf:
        data = {
            "soma_joinid": [0],
        }
        sidf.write(pa.Table.from_pydict(data))

    with tiledbsoma.DataFrame.open(tmp_path.as_posix()) as sidf:
        sidf.read().concat()

    tiledbsoma.stats_dump()
    tiledbsoma.stats_disable()
    stdout, stderr = capfd.readouterr()
    assert stdout != ""
    assert stderr == ""
