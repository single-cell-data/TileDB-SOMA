import pyarrow as pa
import pytest

import tiledbsoma


def test_stats(tmp_path, capsys: pytest.CaptureFixture[str]):
    """Make sure these exist, don't throw, and write correctly."""
    tiledbsoma.tiledbsoma_stats_enable()
    tiledbsoma.tiledbsoma_stats_reset()

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

    tiledbsoma.tiledbsoma_stats_dump()
    tiledbsoma.tiledbsoma_stats_disable()
    stdout, stderr = capsys.readouterr()
    assert stdout != ""
    assert stderr == ""
    print(f"tiledbsoma_stats_dump() = {stdout}")
    tiledbsoma.show_package_versions()


def test_stats_json():
    out = tiledbsoma.tiledbsoma_stats_json()
    assert isinstance(out, str)
    assert out[0] == "["
    assert out[-2] == "]"
    assert out[-1] == "\n"


def test_stats_as_py():
    out = tiledbsoma.tiledbsoma_stats_as_py()
    assert isinstance(out, list)
