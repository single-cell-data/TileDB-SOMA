import tiledbsoma


def test_stats_json():
    out = tiledbsoma.tiledbsoma_stats_json()
    assert isinstance(out, str)
    assert out[0] == "["
    assert out[-2] == "]"
    assert out[-1] == "\n"


def test_stats_as_py():
    out = tiledbsoma.tiledbsoma_stats_as_py()
    assert isinstance(out, list)
