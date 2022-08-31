import os
from pathlib import Path
from typing import Any

import tiledbsc
import tiledbsc.io

HERE = Path(__file__).parent


def soma_obs_callback(soma: tiledbsc.SOMA, _: Any) -> Any:
    output = soma.obs.query(query_string='cell_type == "B cell"', attrs=["cell_type"])
    return list(output.index)


def soma_X_callback(soma: tiledbsc.SOMA, obs_ids: Any) -> Any:
    return soma.X.data.df(obs_ids)


def test_mapper(tmp_path):
    # Populate the collection

    soco_dir = tmp_path.as_posix()

    soco = tiledbsc.SOMACollection(soco_dir)
    soco.create_unless_exists()

    somas = []

    for name, h5ad_path in [
        ("subset-soma-01", HERE.parent / "anndata/subset-soma-01.h5ad"),
        ("subset-soma-02", HERE.parent / "anndata/subset-soma-02.h5ad"),
        ("subset-soma-03", HERE.parent / "anndata/subset-soma-03.h5ad"),
        ("subset-soma-04", HERE.parent / "anndata/subset-soma-04.h5ad"),
    ]:
        soma_path = os.path.join(soco_dir, name)
        soma = tiledbsc.SOMA(soma_path)
        tiledbsc.io.from_h5ad(soma, h5ad_path)
        somas.append(soma)

    for soma in somas:
        soco.add(soma)

    # Do the query
    for parallelize in [False, True]:
        obs_ids_per_soma = soco.map(
            soma_obs_callback,
            None,
            parallelize=parallelize,
        )
        X_dfs = soco.map(
            soma_X_callback,
            obs_ids_per_soma,
            parallelize=parallelize,
        )

        assert len(obs_ids_per_soma) == 4
        assert len(X_dfs) == 4

        assert len(obs_ids_per_soma["subset-soma-01"]) == 7
        assert len(obs_ids_per_soma["subset-soma-02"]) == 7
        assert len(obs_ids_per_soma["subset-soma-03"]) == 4
        assert len(obs_ids_per_soma["subset-soma-04"]) == 6

        assert X_dfs["subset-soma-01"].shape == (15052, 1)
        assert X_dfs["subset-soma-02"].shape == (13391, 1)
        assert X_dfs["subset-soma-03"].shape == (7532, 1)
        assert X_dfs["subset-soma-04"].shape == (8904, 1)
