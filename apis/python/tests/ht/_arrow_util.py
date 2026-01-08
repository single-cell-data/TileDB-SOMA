"""Various utilities for dealing with Arrow data."""

from __future__ import annotations

import pyarrow as pa


def combine_chunks(a: pa.ChunkedArray) -> pa.Array:
    """Semantically identical to pa.ChunkedArray.combine_chunks, but handles the
    `large_` types which are unimplemented by pyarrow.
    """
    type = a.type

    if pa.types.is_large_string(type):
        return a.cast(pa.string()).combine_chunks().cast(type)

    if pa.types.is_large_binary(type):
        return a.cast(pa.binary()).combine_chunks().cast(type)

    if pa.types.is_dictionary(type):
        if pa.types.is_large_string(type.value_type):
            return (
                a
                .cast(
                    pa.dictionary(
                        index_type=type.index_type,
                        value_type=pa.string(),
                        ordered=type.ordered,
                    ),
                )
                .combine_chunks()
                .cast(type)
            )

        if pa.types.is_large_binary(type.value_type):
            return (
                a
                .cast(
                    pa.dictionary(
                        index_type=type.index_type,
                        value_type=pa.binary(),
                        ordered=type.ordered,
                    ),
                )
                .combine_chunks()
                .cast(type)
            )

    return a.combine_chunks()
