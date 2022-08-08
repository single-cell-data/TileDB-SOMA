import os
import pandas as pd
import pyarrow as pa
import libtiledbsc as sc

VERBOSE = True

TEST_DIR = os.path.dirname(__file__)
SOMA_URI = f"{TEST_DIR}/../../test/soco/pbmc3k_processed"

if VERBOSE:
    sc.config_logging("debug")


def test_soma_list():
    soma = sc.SOMA(SOMA_URI)
    array_uris = soma.list_arrays()
    assert len(array_uris) == 19


def test_soma_query():
    soma = sc.SOMA(SOMA_URI)
    sq = soma.query()

    dfs = []
    while chunk := sq.next_results():
        if VERBOSE:
            print("---- CHUNK ----")
            print(chunk)
        dfs.append(chunk["soma/X/data"].to_pandas())

    df = pd.concat(dfs)

    if VERBOSE:
        print("---- FINAL ----")
        print(df)

    assert len(df) == 4848644


def test_soma_buffer_size():
    for i in range(22, 26):
        config = {}
        config["soma.init_buffer_bytes"] = f"{1 << i}"
        sc.debug(f"SOMA query: {config}")

        soma = sc.SOMA(SOMA_URI, config)
        sq = soma.query()

        dfs = []
        while chunk := sq.next_results():
            if VERBOSE:
                print("---- CHUNK ----")
                print(chunk)
            dfs.append(chunk["soma/X/data"].to_pandas())

        df = pd.concat(dfs)

        if VERBOSE:
            print("---- FINAL ----")
            print(df)

        assert len(df) == 4848644


def test_soma_slice_obs():
    soma = sc.SOMA(SOMA_URI)
    sq = soma.query()

    sq.select_obs_attrs(["louvain"])
    sq.set_obs_condition("louvain", "B cells", 4)  # EQ = 4

    dfs = []
    while chunk := sq.next_results():
        if VERBOSE:
            print("---- CHUNK ----")
            print(chunk)
        dfs.append(chunk["soma/X/data"].to_pandas())

    df = pd.concat(dfs)

    if VERBOSE:
        print("---- FINAL ----")
        print(df)

    assert len(df) == 628596


def test_soma_slice_var():
    soma = sc.SOMA(SOMA_URI)
    sq = soma.query()

    sq.select_var_attrs(["n_cells"])
    sq.set_var_condition("n_cells", 50, 0)  # LT = 0

    dfs = []
    while chunk := sq.next_results():
        if VERBOSE:
            print("---- CHUNK ----")
            print(chunk)
        dfs.append(chunk["soma/X/data"].to_pandas())

    df = pd.concat(dfs)

    if VERBOSE:
        print("---- FINAL ----")
        print(df)

    assert len(df) == 1308448


def test_soma_slice_ids():
    soma = sc.SOMA(SOMA_URI)
    sq = soma.query()

    obs_ids = ["AAACATACAACCAC-1", "AAACATTGATCAGC-1", "TTTGCATGCCTCAC-1"]
    var_ids = ["AAGAB", "AAR2", "ZRANB3"]
    sq.select_obs_ids(obs_ids)
    sq.select_var_ids(var_ids)

    dfs = []
    while chunk := sq.next_results():
        if VERBOSE:
            print("---- CHUNK ----")
            print(chunk)
        dfs.append(chunk["soma/X/data"].to_pandas())

    df = pd.concat(dfs)

    if VERBOSE:
        print("---- FINAL ----")
        print(df)

    assert len(df) == 9


if __name__ == "__main__":
    test_soma_query()
