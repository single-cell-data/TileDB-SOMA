#!/usr/bin/env python
import faulthandler
import pathlib
import sys

import pyarrow as pa

import tiledbsoma as soma

faulthandler.enable()


def test():
    path = pathlib.Path("data") / "soma-experiment-versions-2025-04-04" / "1.7.3" / "pbmc3k_processed"

    def read_axis_df(axis_df, coords):
        return axis_df.read(coords).concat(), axis_df.uri

    def read_slot(slot_df, coords):
        return slot_df.read(coords).tables().concat(), slot_df.uri

    with soma.open(path.as_posix()) as exp:
        tp = exp.context.threadpool

        (obs_df, _), (var_df, _) = tp.map(
            read_axis_df, (exp.obs, exp.ms["RNA"].var), ((slice(0, 499),), (slice(None),))
        )
        # obs_joinids = obs_df["soma_joinid"]
        # var_joinids = var_df["soma_joinid"]
        obs_joinids = pa.array(range(500))
        var_joinids = pa.array(range(1838))

        slot_arrays = [
            exp.ms["RNA"][C][K]
            for C, K in [
                ("X", "data"),
                ("obsm", "X_pca"),
                ("obsm", "X_draw_graph_fr"),
                ("obsm", "X_tsne"),
                ("obsm", "X_umap"),
                ("obsp", "connectivities"),
                ("obsp", "distances"),
                ("varm", "PCs"),
            ]
        ]
        futures = [tp.submit(read_slot, slot_df, (obs_joinids, var_joinids)) for slot_df in slot_arrays]
        for ftr in futures:
            data, uri = ftr.result()
            print(f"{uri} complete", file=sys.stderr)
            assert data is not None


def main():
    for _ in range(5000):
        test()
    return 0


if __name__ == "__main__":
    sys.exit(main())
