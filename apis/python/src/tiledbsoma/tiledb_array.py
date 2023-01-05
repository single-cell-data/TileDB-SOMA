from typing import Dict, Optional, Sequence, Union

import tiledb

import tiledbsoma

from .tiledb_object import TileDBObject
from .types import MNTupleStr, MNTupleStrNone


class TileDBArray(TileDBObject):
    """
    Wraps arrays from TileDB-Py by retaining a URI, options, etc.
    Also serves as an abstraction layer to hide TileDB-specific details from the API, unless
    requested.
    """

    def __init__(
        self, uri: str, name: str, *, parent: Optional["tiledbsoma.TileDBGroup"] = None
    ):
        """
        See the TileDBObject constructor.
        """
        super().__init__(uri, name, parent=parent)

    def _open(self, mode: str = "r") -> tiledb.Array:
        """
        This is just a convenience wrapper allowing 'with self._open() as A: ...' rather than
        'with tiledb.open(self.uri) as A: ...'.
        """
        assert mode in ["w", "r"]
        # This works in either 'with self._open() as A:' or 'A = self._open(); ...; A.close().  The
        # reason is that with-as invokes our return value's __enter__ on return from this method,
        # and our return value's __exit__ on exit from the body of the with-block. The tiledb
        # array object does both of those things. (And if it didn't, we'd get a runtime AttributeError
        # on with-as, flagging the non-existence of the __enter__ or __exit__.)
        return tiledb.open(self.uri, mode=mode, ctx=self._ctx)

    def exists(self) -> bool:
        """
        Tells whether or not there is storage for the array. This might be in case a SOMA
        object has not yet been populated, e.g. before calling ``from_anndata`` --- or, if the
        SOMA has been populated but doesn't have this member (e.g. not all SOMAs have a ``varp``).
        """
        with tiledb.scope_ctx(self._ctx):
            return bool(tiledb.array_exists(self.uri))

    def tiledb_array_schema(self) -> tiledb.ArraySchema:
        """
        Returns the TileDB array schema.
        """
        with self._open() as A:
            return A.schema

    def dim_names(self) -> Sequence[str]:
        """
        Reads the dimension names from the schema: for example, ['obs_id', 'var_id'].
        """
        with self._open() as A:
            return [A.schema.domain.dim(i).name for i in range(A.schema.domain.ndim)]

    def dim_names_to_types(self) -> Dict[str, str]:
        """
        Returns a dict mapping from dimension name to dimension type.
        """
        with self._open() as A:
            dom = A.schema.domain
            return {dom.dim(i).name: dom.dim(i).dtype for i in range(dom.ndim)}

    def attr_names(self) -> Sequence[str]:
        """
        Reads the attribute names from the schema: for example, the list of column names in a dataframe.
        """
        with self._open() as A:
            return [A.schema.attr(i).name for i in range(A.schema.nattr)]

    def attr_names_to_types(self) -> Dict[str, str]:
        """
        Returns a dict mapping from attribute name to attribute type.
        """
        with self._open() as A:
            schema = A.schema
            return {
                schema.attr(i).name: schema.attr(i).dtype for i in range(schema.nattr)
            }

    def has_attr_name(self, attr_name: str) -> bool:
        """
        Returns true if the array has the specified attribute name, false otherwise.
        """
        return attr_name in self.attr_names()

    def has_attr_names(self, attr_names: Sequence[str]) -> bool:
        """
        Returns true if the array has all of the specified attribute names, false otherwise.
        """
        attr_names_set = set(self.attr_names())
        return all([attr_name in attr_names_set for attr_name in attr_names])

    def show_metadata(self, recursively: bool = True, indent: str = "") -> None:
        """
        Shows metadata for the array.
        """
        print(f"{indent}[{self.name}]")
        for key, value in self.metadata().items():
            print(f"{indent}- {key}: {value}")

    def _ned_value_to_string(self, value: Union[bytes, str, None]) -> Union[None, str]:
        """
        Helper function for _get_non_empty_domain_as_strings
        """
        if value is None:
            return value
        elif isinstance(value, bytes):
            return value.decode()
        elif isinstance(value, str):
            return value
        else:
            raise NotImplementedError(f"expected bytes or str; got {type(value)}")

    def _get_non_empty_domain_as_strings(self, expected_ndim: int) -> MNTupleStrNone:
        """
        Returns nonempty domain tuples. Due to an implementation detail in TileDB-Py, we get these
        from TileDB-Py as tuples of bytes, not strings. Here we take care of that, as a
        keystroke-saver for callers. This is a helper for resume-ingest mode.
        """
        with self._open() as A:
            # Note that if an array has had its schema written but no data written,
            # its non-empty domain will have `None` values.
            ned = A.nonempty_domain()

        assert len(ned) == expected_ndim
        if expected_ndim == 1:
            lo, hi = ned[0]
            lo = self._ned_value_to_string(lo)
            hi = self._ned_value_to_string(hi)
            return ((lo, hi),)

        elif expected_ndim == 2:
            row_lo, row_hi = ned[0]
            row_lo = self._ned_value_to_string(row_lo)
            row_hi = self._ned_value_to_string(row_hi)
            col_lo, col_hi = ned[1]
            col_lo = self._ned_value_to_string(col_lo)
            col_hi = self._ned_value_to_string(col_hi)
            return ((row_lo, row_hi), (col_lo, col_hi))

        else:
            raise NotImplementedError(
                f"only ndims 1 or 2 are supported; got {expected_ndim}"
            )

    def _chunk_is_contained_in(
        self,
        chunk_mbr: MNTupleStr,
        storage_nonempty_domain: MNTupleStrNone,
    ) -> bool:
        """
        Determines if a dim range is included within the array's non-empty domain.  Note that
        `_get_non_empty_domain_as_string()` is provided to cache non-empty domain values.  Ranges
        are inclusive on both endpoints.  This is a helper for resume-ingest mode.
        """
        assert len(chunk_mbr) == len(storage_nonempty_domain)
        ndim = len(chunk_mbr)

        for i in range(ndim):
            chunk_lo, chunk_hi = chunk_mbr[i]
            storage_lo, storage_hi = storage_nonempty_domain[i]

            if storage_lo is None or storage_hi is None:
                # E.g. an array has had its schema created but no data written yet
                return False

            if chunk_lo < storage_lo or chunk_lo > storage_hi:
                return False
            if chunk_hi < storage_lo or chunk_hi > storage_hi:
                return False

        return True
