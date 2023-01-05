from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Any, Iterator, Mapping, Optional, Sequence, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
import tiledb

import tiledbsoma.util as util

from .logging import log_io, logger
from .tiledb_group import TileDBGroup
from .uns_array import UnsArray


class UnsGroup(TileDBGroup):
    """
    Nominally for soma uns.
    """

    # ----------------------------------------------------------------
    def __init__(self, uri: str, name: str, *, parent: Optional[TileDBGroup] = None):
        """
        See the TileDBObject constructor.
        """
        super().__init__(uri=uri, name=name, parent=parent)

    # ----------------------------------------------------------------
    def keys(self) -> Sequence[str]:
        """
        For uns, ``.keys()`` is a keystroke-saver for the more general group-member
        accessor ``.get_member_names()``.
        """
        return self.get_member_names()

    # At the tiledb-py API level, *all* groups are name-indexable.  But here at the tiledbsoma-py
    # level, we implement name-indexing only for some groups:
    #
    # * Most soma member references are done using Python's dot syntax. For example, rather than
    #   soma['X'], we have simply soma.X, and likewise, soma.raw.X.  Likewise soma.obs and soma.var.
    #
    # * Index references are supported for obsm, varm, obsp, varp, and uns. E.g.
    #   soma.obsm['X_pca'] or soma.uns['neighbors']['params']['method']
    #
    # * Overloading the ``[]`` operator at the TileDBGroup level isn't necessary -- e.g. we don't need
    #   soma['X'] when we have soma.X -- but also it causes circular-import issues in Python.
    #
    # * Rather than doing a TileDBIndexableGroup which overloads the ``[]`` operator, we overload
    #   the ``[]`` operator separately in the various classes which need indexing. This is again to
    #   avoid circular-import issues, and means that [] on ``AnnotationMatrixGroup`` will return an
    #   ``AnnotationMatrix, [] on ``UnsGroup`` will return ``UnsArray`` or ``UnsGroup``, etc.
    def __getitem__(self, name: str) -> Union[UnsGroup, UnsArray, None]:
        """
        Returns an ``UnsArray`` or ``UnsGroup`` element at the given name within the group, or None if
        no such member exists.  Overloads the [...] operator.
        """

        with self._open("r") as G:
            if name not in G:
                return None

            obj = G[name]  # This returns a tiledb.object.Object.
            if obj.type == tiledb.tiledb.Group:
                return UnsGroup(uri=obj.uri, name=name, parent=self)
            elif obj.type == tiledb.libtiledb.Array:
                return UnsArray(uri=obj.uri, name=name, parent=self)
            else:
                raise Exception(
                    f"Internal error: found uns group element neither subgroup nor array: type is {str(obj.type)}"
                )

    # ----------------------------------------------------------------
    def __contains__(self, name: str) -> bool:
        """
        Implements '"namegoeshere" in soma.uns'.
        """
        with self._open("r") as G:
            return name in G

    # ----------------------------------------------------------------
    def __iter__(self) -> Iterator[Union[UnsGroup, UnsArray]]:
        """
        Implements ``for element in soma.uns: ...``
        """
        with self._open("r") as G:
            for obj in G:  # tiledb.object.Object
                if obj.type == tiledb.tiledb.Group:
                    yield UnsGroup(uri=obj.uri, name=obj.name, parent=self)
                elif obj.type == tiledb.libtiledb.Array:
                    yield UnsArray(uri=obj.uri, name=obj.name, parent=self)
                else:
                    raise Exception(
                        f"Internal error: found uns group element neither subgroup nor array: type is {obj.type}"
                    )

    # ----------------------------------------------------------------
    def __repr__(self) -> str:
        """
        Default display for uns groups.
        """
        return self._repr_aux()

    # ----------------------------------------------------------------
    def _repr_aux(self, display_name: str = "uns", indent: str = "") -> str:
        """
        Recursively displays the uns data.
        """
        strings = []
        strings.append(indent + display_name + ":")
        # Scalars are stored as metadata, not 1D arrays
        for k, v in self.metadata().items():
            if not k.startswith("__"):
                strings.append(indent + k + ": " + repr(v))
        # Now do subgroups, and non-scalar values
        for e in self:
            element_display_name = display_name + "/" + e.name
            if isinstance(e, UnsGroup):
                strings.append(
                    indent
                    + e._repr_aux(
                        display_name=element_display_name, indent=indent + "  "
                    )
                )
            else:
                strings.append(indent + element_display_name + "/")
                with e._open() as A:
                    strings.append(indent + str(A[:]))
        return "\n".join(strings)

    # ----------------------------------------------------------------
    def _from_anndata_uns_aux(
        self, key: str, value: Any, component_uri: str, ingest_mode: str
    ) -> None:
        """
        Helper method for `from_anndata_uns`.
        """

        if key == "rank_genes_groups":
            # TODO:
            # This is of type 'structured array':
            # https://numpy.org/doc/stable/user/basics.rec.html
            #
            # >>> a.uns['rank_genes_groups']['names'].dtype
            # dtype([('0', 'O'), ('1', 'O'), ('2', 'O'), ('3', 'O'), ('4', 'O'), ('5', 'O'), ('6', 'O'), ('7', 'O')])
            # >>> type(a.uns['rank_genes_groups']['names'])
            # <class 'numpy.ndarray'>
            #
            # We don’t have a way to model this directly in TileDB schema right now. We support
            # multiplicities of a single scalar type, e.g. a record array with cell_val_num==3
            # and float32 slots (which would correspond to numpy record array
            # np.dtype([("field1", "f4"), ("field2", "f4"), ("field3", "f4",)])). We don’t
            # support nested cells, AKA "list" type.
            #
            # This could, however, be converted to a dataframe and ingested that way.
            log_io(None, f"{self._indent}Skipping structured array: {component_uri}")
            return

        if isinstance(value, Mapping):
            # Nested data, e.g. a.uns['draw-graph']['params']['layout']
            subgroup = UnsGroup(uri=component_uri, name=key, parent=self)
            subgroup.from_anndata_uns(value, ingest_mode)
            self._add_object(subgroup)
            return

        # Write scalars as metadata, not length-1 component arrays
        if isinstance(value, np.str_):
            with self._open("w") as G:
                G.meta[key] = value
            # TODO: WUT
            # Needs explicit cast from numpy.str_ to str for tiledb.from_numpy
            return

        if isinstance(value, (int, float, str)):
            # Nominally this is unit-test data
            with self._open("w") as G:
                G.meta[key] = value
            return

        # Everything else is a component array, or unhandleable
        array = UnsArray(uri=component_uri, name=key, parent=self)

        if ingest_mode == "resume" and array.exists():
            log_io(
                f"Skipped {array.nested_name}",
                f"{self._indent}Skipping existing array {array.nested_name}",
            )
            return

        if isinstance(value, pd.DataFrame):
            array.from_pandas_dataframe(value)
            self._add_object(array)
            return

        if isinstance(value, sp.csr_matrix):
            array.from_scipy_csr(value)
            self._add_object(array)
            return

        if array._maybe_from_numpyable_object(value):
            self._add_object(array)
            return

        logger.error(
            f"{self._indent}Skipping unrecognized type: {component_uri} {type(value)}",
        )

    def from_anndata_uns(self, uns: Mapping[str, Any], ingest_mode: str) -> None:
        """
        Populates the uns group for the soma object.

        :param uns: anndata.uns.
        """

        s = util.get_start_stamp()
        log_io(None, f"{self._indent}START  WRITING {self.nested_name}")

        # Must be done first, to create the parent directory
        self.create_unless_exists()

        # See comments in _get_child_uris
        child_uris = self._get_child_uris(list(uns.keys()))

        futures = []
        max_workers = self._soma_options.max_thread_pool_workers

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            for key in uns.keys():
                future = executor.submit(
                    self._from_anndata_uns_aux,
                    key,
                    uns[key],
                    child_uris[key],
                    ingest_mode,
                )
                futures.append(future)

        for future in futures:
            future.result()

        log_io(
            None,  # At info level, print only arrays not groups
            util.format_elapsed(s, f"{self._indent}FINISH WRITING {self.nested_name}"),
        )

    # ----------------------------------------------------------------
    def to_dict_of_matrices(self) -> Mapping[str, Any]:
        """
        Reads the recursive group/array uns data from TileDB storage
        and returns them as a recursive dict of matrices.
        """
        if not self.exists():
            log_io(
                f"{self._indent}{self.nested_name} not found",
                f"{self._indent}{self.nested_name} not found",
            )
            return {}

        s = util.get_start_stamp()
        log_io(None, f"{self._indent}START  read {self.nested_name}")

        with self._open() as G:
            retval = {}
            for element in G:
                name = element.name

                if element.type == tiledb.tiledb.Group:
                    child_group = UnsGroup(uri=element.uri, name=name, parent=self)
                    retval[name] = child_group.to_dict_of_matrices()

                elif element.type == tiledb.libtiledb.Array:
                    child_array = UnsArray(uri=element.uri, name=name, parent=self)
                    retval[name] = child_array.to_matrix()

                else:
                    raise Exception(
                        f"Internal error: found uns group element neither group nor array: type is {str(element.type)}"
                    )

        log_io(
            f"Read {self.nested_name}",
            util.format_elapsed(s, f"{self._indent}FINISH READING {self.nested_name}"),
        )

        return retval
