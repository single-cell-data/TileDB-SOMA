import os
import pandas as pd
import pyarrow as pa
import libtiledbsc as sc

VERBOSE = False

TEST_DIR = os.path.dirname(__file__)
SOMA_URI = f"{TEST_DIR}/../../test/soco/pbmc3k_processed"

if VERBOSE:
    sc.config_logging("debug")


def test_soma_list():
    soma = sc.SOMA(SOMA_URI)
    array_uris = soma.list_arrays()
    assert len(array_uris) == 19


# TODO: Debug this message from pa.concat_tables
# var_id: [<Invalid array: Length spanned by binary offsets (6010916) larger than values array (size 6010627)>,
def test_soma_query():
    soma = sc.SOMA(SOMA_URI)
    sq = soma.query()

    table = sq.next_results()
    while chunk := sq.next_results():
        if VERBOSE:
            print("---- CHUNK ----")
            print(chunk)
        table = pa.concat_tables([table, chunk])

    if VERBOSE:
        print("---- FINAL ----")
        print(table)

    assert len(table) == 4848644


def test_soma_buffer_size():
    for i in range(22, 26):
        config = {}
        config["soma.init_buffer_bytes"] = f"{1 << i}"
        sc.debug(f"SOMA query: {config}")
        soma = sc.SOMA(SOMA_URI, config)
        sq = soma.query()

        table = sq.next_results()
        while chunk := sq.next_results():
            table = pa.concat_tables([table, chunk])

        if VERBOSE:
            print(table)

        assert len(table) == 4848644


def test_soma_slice_obs():
    soma = sc.SOMA(SOMA_URI)
    sq = soma.query()

    sq.select_obs_attrs(["louvain"])
    sq.set_obs_condition("louvain", "B cells", 4)  # EQ = 4

    table = sq.next_results()
    while chunk := sq.next_results():
        table = pa.concat_tables([table, chunk])

    assert len(table) == 628596


def test_soma_slice_var():
    config = {}
    # config["config.logging_level"] = "5"
    # config["soma.init_buffer_bytes"] = f"{1 << 26}"

    soma = sc.SOMA(SOMA_URI, config)
    sq = soma.query()

    sq.select_var_attrs(["n_cells"])
    sq.set_var_condition("n_cells", 50, 0)  # LT = 0

    table = sq.next_results()
    while chunk := sq.next_results():
        table = pa.concat_tables([table, chunk])

    assert len(table) == 1308448


def test_soma_slice_ids():
    soma = sc.SOMA(SOMA_URI)
    sq = soma.query()

    obs_ids = ["AAACATACAACCAC-1", "AAACATTGATCAGC-1", "TTTGCATGCCTCAC-1"]
    var_ids = ["AAGAB", "AAR2", "ZRANB3"]
    sq.select_obs_ids(obs_ids)
    sq.select_var_ids(var_ids)

    table = sq.next_results()
    while chunk := sq.next_results():
        table = pa.concat_tables([table, chunk])

    assert len(table) == 9


if __name__ == "__main__":
    test_soma_query()
