from typing import Dict, List

import tiledb

import tiledbsc.util_tiledb

from .tiledb_object import TileDBObject


class TileDBArray(TileDBObject):
    """
    Wraps arrays from TileDB-Py by retaining a URI, options, etc.
    Also serves as an abstraction layer to hide TileDB-specific details from the API, unless
    requested.
    """

    def __init__(
        self,
        uri: str,
        name: str,
        # It's a circular import if we say this, but this is really:
        # parent: Optional[TileDBGroup] = None,
        parent=None,
    ):
        """
        See the TileDBObject constructor.
        """
        super().__init__(uri=uri, name=name, parent=parent)

    def _object_type(self) -> str:
        """
        This should be implemented by child classes and should return what tiledb.object_type(uri)
        returns for objects of a given type -- nominally 'group' or 'array'.
        """
        return "array"

    def _open(self, mode="r"):
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
        object has not yet been populated, e.g. before calling `from_anndata` -- or, if the
        SOMA has been populated but doesn't have this member (e.g. not all SOMAs have a `varp`).
        """
        return tiledb.array_exists(self.uri)

    def tiledb_array_schema(self):
        """
        Returns the TileDB array schema.
        """
        with self._open() as A:
            return A.schema

    def dim_names(self) -> List[str]:
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

    def attr_names(self) -> List[str]:
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

    def has_attr_names(self, attr_names: List[str]) -> bool:
        """
        Returns true if the array has all of the specified attribute names, false otherwise.
        """
        attr_names_set = set(self.attr_names())
        return all([attr_name in attr_names_set for attr_name in attr_names])

    def _set_object_type_metadata(self) -> None:
        """
        This helps nested-structured traversals (especially those that start at the SOMACollection
        level) confidently navigate with a minimum of introspection on group contents.
        """
        with self._open("w") as A:
            A.meta[
                tiledbsc.util.SOMA_OBJECT_TYPE_METADATA_KEY
            ] = self.__class__.__name__
            A.meta[
                tiledbsc.util.SOMA_ENCODING_VERSION_METADATA_KEY
            ] = tiledbsc.util.SOMA_ENCODING_VERSION

    def get_object_type(self) -> str:
        """
        Returns the class name associated with the array.
        """
        with self._open("r") as A:
            return A.meta[tiledbsc.util_tiledb.SOMA_OBJECT_TYPE_METADATA_KEY]

    def show_metadata(self, recursively=True, indent="") -> None:
        """
        Shows metadata for the array.
        """
        print(f"{indent}[{self.name}]")
        for key, value in self.metadata().items():
            print(f"{indent}- {key}: {value}")
