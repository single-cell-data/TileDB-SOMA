import pandas as pd
import pyarrow as pa
import libtiledbsc as sc

SOMA_URI = "/home/gspowley/work/TileDB-SingleCell/test/soco/pbmc3k_processed"


def test_soma_list():
    soma = sc.SOMA(SOMA_URI)
    array_uris = soma.list_arrays()
    assert len(array_uris) == 19


def test_soma_query():
    for i in range(22, 26):
        config = {"soma.init_buffer_bytes": f"{1 << i}"}
        soma = sc.SOMA(SOMA_URI, config)
        sq = soma.query()

        table = sq.next_results()
        while chunk := sq.next_results():
            table = pa.concat_tables([table, chunk])

        assert len(table.to_pandas()) == 4848644


if __name__ == "__main__":
    test_soma()
