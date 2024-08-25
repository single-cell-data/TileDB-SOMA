# Copyright (c) 2021-2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2024 TileDB, Inc.
#
# Licensed under the MIT License.

from __future__ import annotations

import pathlib
from typing import Callable, List, Optional, Sequence, Union
from unittest.mock import patch

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from scipy import sparse
from scipy.sparse import coo_matrix, spmatrix

import tiledbsoma as soma
from tiledbsoma import Experiment, _factory
from tiledbsoma._collection import CollectionBase

# conditionally import torch, as it will not be available in all test environments.
# This supports the pytest `ml` mark, which can be used to disable all PyTorch-dependent
# tests.
try:
    from tiledbsoma_ml.pytorch import (
        ExperimentAxisQueryDataPipe,
        ExperimentAxisQueryIterableDataset,
        ExperimentAxisQueryIterable,
        experiment_dataloader,
    )
    from torch.utils.data._utils.worker import WorkerInfo
except ImportError:
    # this should only occur when not running `ml`-marked tests
    pass


def pytorch_x_value_gen(obs_range: range, var_range: range) -> spmatrix:
    occupied_shape = (
        obs_range.stop - obs_range.start,
        var_range.stop - var_range.start,
    )
    checkerboard_of_ones = coo_matrix(np.indices(occupied_shape).sum(axis=0) % 2)
    checkerboard_of_ones.row += obs_range.start
    checkerboard_of_ones.col += var_range.start
    return checkerboard_of_ones


def pytorch_seq_x_value_gen(obs_range: range, var_range: range) -> spmatrix:
    """A sparse matrix where the values of each col are the obs_range values. Useful for checking the
    X values are being returned in the correct order."""
    data = np.vstack([list(obs_range)] * len(var_range)).flatten()
    rows = np.vstack([list(obs_range)] * len(var_range)).flatten()
    cols = np.column_stack([list(var_range)] * len(obs_range)).flatten()
    return coo_matrix((data, (rows, cols)))


@pytest.fixture
def X_layer_names() -> List[str]:
    return ["raw"]


@pytest.fixture
def obsp_layer_names() -> Optional[List[str]]:
    return None


@pytest.fixture
def varp_layer_names() -> Optional[List[str]]:
    return None


def add_dataframe(coll: CollectionBase, key: str, value_range: range) -> None:
    df = coll.add_new_dataframe(
        key,
        schema=pa.schema(
            [
                ("soma_joinid", pa.int64()),
                ("label", pa.large_string()),
                ("label2", pa.large_string()),
            ]
        ),
        index_column_names=["soma_joinid"],
    )
    df.write(
        pa.Table.from_pydict(
            {
                "soma_joinid": list(value_range),
                "label": [str(i) for i in value_range],
                "label2": ["c" for i in value_range],
            }
        )
    )


def add_sparse_array(
    coll: CollectionBase,
    key: str,
    obs_range: range,
    var_range: range,
    value_gen: Callable[[range, range], spmatrix],
) -> None:
    a = coll.add_new_sparse_ndarray(
        key, type=pa.float32(), shape=(obs_range.stop, var_range.stop)
    )
    tensor = pa.SparseCOOTensor.from_scipy(value_gen(obs_range, var_range))
    a.write(tensor)


@pytest.fixture(scope="function")
def soma_experiment(
    tmp_path: pathlib.Path,
    obs_range: Union[int, range],
    var_range: Union[int, range],
    X_value_gen: Callable[[range, range], sparse.spmatrix],
    obsp_layer_names: Sequence[str],
    varp_layer_names: Sequence[str],
) -> soma.Experiment:
    with soma.Experiment.create((tmp_path / "exp").as_posix()) as exp:
        if isinstance(obs_range, int):
            obs_range = range(obs_range)
        if isinstance(var_range, int):
            var_range = range(var_range)

        add_dataframe(exp, "obs", obs_range)
        ms = exp.add_new_collection("ms")
        rna = ms.add_new_collection("RNA", soma.Measurement)
        add_dataframe(rna, "var", var_range)
        rna_x = rna.add_new_collection("X", soma.Collection)
        add_sparse_array(rna_x, "raw", obs_range, var_range, X_value_gen)

        if obsp_layer_names:
            obsp = rna.add_new_collection("obsp")
            for obsp_layer_name in obsp_layer_names:
                add_sparse_array(
                    obsp, obsp_layer_name, obs_range, var_range, X_value_gen
                )

        if varp_layer_names:
            varp = rna.add_new_collection("varp")
            for varp_layer_name in varp_layer_names:
                add_sparse_array(
                    varp, varp_layer_name, obs_range, var_range, X_value_gen
                )
    return _factory.open((tmp_path / "exp").as_posix())


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(6, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test_non_batched(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
    use_eager_fetch: bool,
) -> None:
    # batch_size should default to 1
    exp_data_pipe = PipeClass(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    assert type(exp_data_pipe.shape) is tuple
    assert len(exp_data_pipe.shape) == 2
    assert exp_data_pipe.shape == (6, 3)

    row_iter = iter(exp_data_pipe)

    row = next(row_iter)
    assert isinstance(row[0], np.ndarray)
    assert isinstance(row[1], pd.DataFrame)
    assert row[0].shape == (3,)
    assert row[1].shape == (1, 1)
    assert row[0].tolist() == [0, 1, 0]
    assert row[1].keys() == ["label"]
    assert row[1]["label"].tolist() == ["0"]


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(6, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test_uneven_soma_and_result_batches(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
    use_eager_fetch: bool,
) -> None:
    """This is checking that batches are correctly created when they require fetching multiple chunks."""
    exp_data_pipe = PipeClass(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        shuffle=False,
        batch_size=3,
        io_batch_size=2,
        use_eager_fetch=use_eager_fetch,
    )
    row_iter = iter(exp_data_pipe)

    X_batch, obs_batch = next(row_iter)
    assert isinstance(X_batch, np.ndarray)
    assert isinstance(obs_batch, pd.DataFrame)
    assert X_batch.shape[0] == obs_batch.shape[0]
    assert X_batch.shape == (3, 3)
    assert obs_batch.shape == (3, 1)
    assert X_batch[0].tolist() == [0, 1, 0]
    assert ["label"] == obs_batch.keys()
    assert obs_batch["label"].tolist() == ["0", "1", "2"]


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(6, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test_batching__all_batches_full_size(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
    use_eager_fetch: bool,
) -> None:
    exp_data_pipe = PipeClass(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        batch_size=3,
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    batch_iter = iter(exp_data_pipe)

    batch = next(batch_iter)
    assert batch[0].tolist() == [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
    assert batch[1].keys() == ["label"]
    assert batch[1]["label"].tolist() == ["0", "1", "2"]

    batch = next(batch_iter)
    assert batch[0].tolist() == [[1, 0, 1], [0, 1, 0], [1, 0, 1]]
    assert batch[1].keys() == ["label"]
    assert batch[1]["label"].tolist() == ["3", "4", "5"]

    with pytest.raises(StopIteration):
        next(batch_iter)


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [
        (range(100_000_000, 100_000_003), 3, pytorch_x_value_gen, use_eager_fetch)
        for use_eager_fetch in (True, False)
    ],
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test_unique_soma_joinids(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
    use_eager_fetch: bool,
) -> None:
    exp_data_pipe = PipeClass(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["soma_joinid", "label"],
        batch_size=3,
        use_eager_fetch=use_eager_fetch,
    )

    soma_joinids = np.concatenate(
        [batch[1]["soma_joinid"].to_numpy() for batch in exp_data_pipe]
    )
    assert len(np.unique(soma_joinids)) == len(soma_joinids)


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(5, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test_batching__partial_final_batch_size(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
    use_eager_fetch: bool,
) -> None:
    exp_data_pipe = PipeClass(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        batch_size=3,
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    batch_iter = iter(exp_data_pipe)

    next(batch_iter)
    batch = next(batch_iter)
    assert batch[0].tolist() == [[1, 0, 1], [0, 1, 0]]

    with pytest.raises(StopIteration):
        next(batch_iter)


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(3, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test_batching__exactly_one_batch(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
    use_eager_fetch: bool,
) -> None:
    exp_data_pipe = PipeClass(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        batch_size=3,
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    batch_iter = iter(exp_data_pipe)

    batch = next(batch_iter)
    assert batch[0].tolist() == [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
    assert batch[1]["label"].tolist() == ["0", "1", "2"]

    with pytest.raises(StopIteration):
        next(batch_iter)


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(6, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test_batching__empty_query_result(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
    use_eager_fetch: bool,
) -> None:
    exp_data_pipe = PipeClass(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_query=soma.AxisQuery(coords=([],)),
        obs_column_names=["label"],
        batch_size=3,
        use_eager_fetch=use_eager_fetch,
    )
    batch_iter = iter(exp_data_pipe)

    with pytest.raises(StopIteration):
        next(batch_iter)


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [
        (10, 1, pytorch_x_value_gen, use_eager_fetch)
        for use_eager_fetch in (True, False)
    ],
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test_batching__partial_soma_batches_are_concatenated(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
    use_eager_fetch: bool,
) -> None:
    exp_data_pipe = PipeClass(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        batch_size=3,
        # set SOMA batch read size such that PyTorch batches will span the tail and head of two SOMA batches
        io_batch_size=4,
        use_eager_fetch=use_eager_fetch,
    )

    full_result = list(exp_data_pipe)

    assert [len(batch[0]) for batch in full_result] == [3, 3, 3, 1]


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen", [(6, 3, pytorch_x_value_gen)]
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test_multiprocessing__returns_full_result(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
) -> None:
    """Tests the ExperimentAxisQueryDataPipe provides all data, as collected from multiple processes that are managed by a
    PyTorch DataLoader with multiple workers configured."""

    dp = PipeClass(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["soma_joinid", "label"],
        io_batch_size=3,  # two chunks, one per worker
    )
    # Note we're testing the ExperimentAxisQueryDataPipe via a DataLoader, since this is what sets up the multiprocessing
    dl = experiment_dataloader(dp, num_workers=2)

    full_result = list(iter(dl))

    soma_joinids = np.concatenate([t[1]["soma_joinid"].to_numpy() for t in full_result])
    assert sorted(soma_joinids) == list(range(6))


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen", [(6, 3, pytorch_x_value_gen)]
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test_distributed__returns_data_partition_for_rank(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
) -> None:
    """Tests pytorch._partition_obs_joinids() behavior in a simulated PyTorch distributed processing mode,
    using mocks to avoid having to do real PyTorch distributed setup."""

    with patch("torch.distributed.is_initialized") as mock_dist_is_initialized, patch(
        "torch.distributed.get_rank"
    ) as mock_dist_get_rank, patch(
        "torch.distributed.get_world_size"
    ) as mock_dist_get_world_size:
        mock_dist_is_initialized.return_value = True
        mock_dist_get_rank.return_value = 1
        mock_dist_get_world_size.return_value = 3

        dp = PipeClass(
            soma_experiment,
            measurement_name="RNA",
            X_name="raw",
            obs_column_names=["soma_joinid"],
            io_batch_size=2,
            shuffle=False,
        )
        full_result = list(iter(dp))

        soma_joinids = np.concatenate(
            [t[1]["soma_joinid"].to_numpy() for t in full_result]
        )

        # Of the 6 obs rows, the PyTorch process of rank 1 should get [2, 3]
        # (rank 0 gets [0, 1], rank 2 gets [4, 5])
        assert sorted(soma_joinids) == [2, 3]


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen", [(12, 3, pytorch_x_value_gen)]
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test_distributed_and_multiprocessing__returns_data_partition_for_rank(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
) -> None:
    """Tests pytorch._partition_obs_joinids() behavior in a simulated PyTorch distributed processing mode and
    DataLoader multiprocessing mode, using mocks to avoid having to do distributed pytorch
    setup or real DataLoader multiprocessing."""

    with patch("torch.utils.data.get_worker_info") as mock_get_worker_info, patch(
        "torch.distributed.is_initialized"
    ) as mock_dist_is_initialized, patch(
        "torch.distributed.get_rank"
    ) as mock_dist_get_rank, patch(
        "torch.distributed.get_world_size"
    ) as mock_dist_get_world_size:
        mock_get_worker_info.return_value = WorkerInfo(id=1, num_workers=2, seed=1234)
        mock_dist_is_initialized.return_value = True
        mock_dist_get_rank.return_value = 1
        mock_dist_get_world_size.return_value = 3

        dp = PipeClass(
            soma_experiment,
            measurement_name="RNA",
            X_name="raw",
            obs_column_names=["soma_joinid"],
            io_batch_size=2,
            shuffle=False,
        )

        full_result = list(iter(dp))

        soma_joinids = np.concatenate(
            [t[1]["soma_joinid"].to_numpy() for t in full_result]
        )

        # Of the 12 obs rows, the PyTorch process of rank 1 should get [4..7], and then within that partition,
        # the 2nd DataLoader process should get the second half of the rank's partition, which is just [6, 7]
        # (rank 0 gets [0..3], rank 2 gets [8..11])
        assert sorted(soma_joinids) == [6, 7]


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(3, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test_experiment_dataloader__non_batched(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
    use_eager_fetch: bool,
) -> None:
    dp = PipeClass(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    dl = experiment_dataloader(dp)
    data = [row for row in dl]
    assert all(d[0].shape == (3,) for d in data)
    assert all(d[1].shape == (1, 1) for d in data)

    row = data[0]
    assert row[0].tolist() == [0, 1, 0]
    assert row[1]["label"].tolist() == ["0"]


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(6, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test_experiment_dataloader__batched(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
    use_eager_fetch: bool,
) -> None:
    dp = PipeClass(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        batch_size=3,
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    dl = experiment_dataloader(dp)
    data = [row for row in dl]

    batch = data[0]
    assert batch[0].tolist() == [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
    assert batch[1].to_numpy().tolist() == [[0], [1], [2]]


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [
        (10, 3, pytorch_x_value_gen, use_eager_fetch)
        for use_eager_fetch in (True, False)
    ],
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test_experiment_dataloader__batched_length(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
    use_eager_fetch: bool,
) -> None:
    dp = PipeClass(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        batch_size=3,
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    dl = experiment_dataloader(dp)
    assert len(dl) == len(list(dl))


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(6, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test__X_tensor_dtype_matches_X_matrix(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
    use_eager_fetch: bool,
) -> None:
    dp = PipeClass(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        batch_size=3,
        use_eager_fetch=use_eager_fetch,
    )
    data = next(iter(dp))

    assert data[0].dtype == np.float32


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen", [(10, 1, pytorch_x_value_gen)]
)
def test__pytorch_splitting(
    soma_experiment: Experiment,
) -> None:
    dp = ExperimentAxisQueryDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
    )
    # function not available for IterableDataset, yet....
    dp_train, dp_test = dp.random_split(weights={"train": 0.7, "test": 0.3}, seed=1234)
    dl = experiment_dataloader(dp_train)

    all_rows = list(iter(dl))
    assert len(all_rows) == 7


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen", [(16, 1, pytorch_seq_x_value_gen)]
)
@pytest.mark.parametrize(
    "PipeClass", (ExperimentAxisQueryDataPipe, ExperimentAxisQueryIterableDataset)
)
def test__shuffle(
    PipeClass: ExperimentAxisQueryDataPipe | ExperimentAxisQueryIterableDataset,
    soma_experiment: Experiment,
) -> None:
    dp = PipeClass(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        shuffle=True,
    )

    all_rows = list(iter(dp))
    assert all(r[0].shape == (1,) for r in all_rows)
    soma_joinids = [row[1]["soma_joinid"].iloc[0] for row in all_rows]
    X_values = [row[0][0].item() for row in all_rows]

    # same elements
    assert set(soma_joinids) == set(range(16))
    # not ordered! (...with a `1/16!` probability of being ordered)
    assert soma_joinids != list(range(16))
    # randomizes X in same order as obs
    # note: X values were explicitly set to match obs_joinids to allow for this simple assertion
    assert X_values == soma_joinids


@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen", [(6, 3, pytorch_x_value_gen)]
)
def test_experiment_axis_query_iterable_error_checks(soma_experiment: Experiment) -> None:
    dp = ExperimentAxisQueryIterable(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        shuffle=True,
    )
    with pytest.raises(NotImplementedError):
        dp[0]

    with pytest.raises(ValueError):
        dp = ExperimentAxisQueryIterable(
            soma_experiment,
            obs_column_names=(),
            measurement_name="RNA",
            X_name="raw",
            shuffle=True,
        )



def test_experiment_dataloader__unsupported_params__fails() -> None:
    with patch(
        "tiledbsoma_ml.pytorch.ExperimentAxisQueryDataPipe"
    ) as dummy_exp_data_pipe:
        with pytest.raises(ValueError):
            experiment_dataloader(dummy_exp_data_pipe, shuffle=True)
        with pytest.raises(ValueError):
            experiment_dataloader(dummy_exp_data_pipe, batch_size=3)
        with pytest.raises(ValueError):
            experiment_dataloader(dummy_exp_data_pipe, batch_sampler=[])
        with pytest.raises(ValueError):
            experiment_dataloader(dummy_exp_data_pipe, sampler=[])
        with pytest.raises(ValueError):
            experiment_dataloader(dummy_exp_data_pipe, collate_fn=lambda x: x)


def test_batched() -> None:
    from tiledbsoma_ml.pytorch import _batched

    assert list(_batched(range(6), 1)) == list((i,) for i in range(6))
    assert list(_batched(range(6), 2)) == [(0, 1), (2, 3), (4, 5)]
    assert list(_batched(range(6), 3)) == [(0, 1, 2), (3, 4, 5)]
    assert list(_batched(range(6), 4)) == [(0, 1, 2, 3), (4, 5)]
    assert list(_batched(range(6), 5)) == [(0, 1, 2, 3, 4), (5,)]
    assert list(_batched(range(6), 6)) == [(0, 1, 2, 3, 4, 5)]
    assert list(_batched(range(6), 7)) == [(0, 1, 2, 3, 4, 5)]

    # bogus batch value
    with pytest.raises(ValueError):
        list(_batched([0, 1], 0))
    with pytest.raises(ValueError):
        list(_batched([2, 3], -1))


def test_splits() -> None:
    from tiledbsoma_ml.pytorch import _splits

    assert _splits(10, 1).tolist() == [0, 10]
    assert _splits(10, 3).tolist() == [0, 4, 7, 10]
    assert _splits(10, 4).tolist() == [0, 3, 6, 8, 10]
    assert _splits(10, 10).tolist() == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    assert _splits(10, 11).tolist() == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10]

    # bad number of sections
    with pytest.raises(ValueError):
        _splits(10, 0)
    with pytest.raises(ValueError):
        _splits(10, -1)


def test_csr_to_dense() -> None:
    from tiledbsoma_ml.pytorch import _csr_to_dense

    coo = sparse.eye(1001, 77, format="coo", dtype=np.float32)

    assert np.array_equal(
        sparse.csr_array(coo).todense(), _csr_to_dense(sparse.csr_array(coo))
    )
    assert np.array_equal(
        sparse.csr_matrix(coo).todense(), _csr_to_dense(sparse.csr_matrix(coo))
    )

    csr = sparse.csr_array(coo)
    assert np.array_equal(csr.todense(), _csr_to_dense(csr))
    assert np.array_equal(csr[1:, :].todense(), _csr_to_dense(csr[1:, :]))
    assert np.array_equal(csr[:, 1:].todense(), _csr_to_dense(csr[:, 1:]))
    assert np.array_equal(csr[3:501, 1:22].todense(), _csr_to_dense(csr[3:501, 1:22]))
