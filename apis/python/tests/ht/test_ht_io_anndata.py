from __future__ import annotations

import pathlib

import anndata
import hypothesis as ht
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra import numpy as ht_np
from hypothesis.extra import pandas as ht_pd

import tiledbsoma as soma
import tiledbsoma.io

from tests.ht._ht_util import posix_filename

# from_anndata(experiment_uri: 'str', anndata: 'ad.AnnData', measurement_name: 'str', *, context: 'SOMATileDBContext | None' = None, platform_config: 'PlatformConfig | None' = None, obs_id_name: 'str' = 'obs_id', var_id_name: 'str' = 'var_id', X_layer_name: 'str' = 'data', raw_X_layer_name: 'str' = 'data', ingest_mode: 'IngestMode' = 'write', use_relative_uri: 'bool | None' = None, X_kind: 'Union[Type[SparseNDArray], Type[DenseNDArray]]' = <class 'tiledbsoma._sparse_nd_array.SparseNDArray'>, registration_mapping: 'ExperimentAmbientLabelMapping | None' = None, uns_keys: 'Sequence[str] | None' = None, additional_metadata: 'AdditionalMetadata' = None) -> 'str'


@st.composite
def anndatas(draw: st.DrawFn) -> anndata.AnnData:
    """
    is empty OK?
    is empty obs/var OK?
    etc
    """

    obs = draw(
        ht_pd.data_frames(columns=columns, index=ht_pd.range_indexes(min_size=min_size))
    )

    return anndata.AnnData()


@settings(suppress_health_check=(ht.HealthCheck.function_scoped_fixture,))
@given(
    data=st.data(),
    adata=st.from_type(anndata.AnnData),
    measurement_name=st.text(min_size=1),
    X_layer_name=st.text(min_size=1),
    raw_X_layer_name=st.text(min_size=1),
)
def test_roundtrip_from_anndata_to_anndata(
    data: st.DataFn,
    adata: anndata.AnnData,
    measurement_name: str,
    X_layer_name: str,
    raw_X_layer_name: str,
    tmp_path_factory,  # fixure
) -> None:

    experiment_uri = tmp_path_factory.mktemp("anndata-").as_posix()

    print("from_anndata(")
    print(adata)
    print(experiment_uri)
    print(measurement_name)
    print(X_layer_name)
    print(raw_X_layer_name)
    print(")")
