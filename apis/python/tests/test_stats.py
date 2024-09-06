import tiledbsoma


def test_stats_json():
    out = tiledbsoma.tiledbsoma_stats_json()
    assert isinstance(out, str)
    assert out[0] == "["
    assert out[-2] == "]"
    assert out[-1] == "\n"


def test_stats_parsed():
    out = tiledbsoma.tiledbsoma_stats_parsed()
    assert isinstance(out, list)
