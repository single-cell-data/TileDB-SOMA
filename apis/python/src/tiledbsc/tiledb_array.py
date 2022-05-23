import tiledb
from .soma_options import SOMAOptions
from .tiledb_object import TileDBObject
from .tiledb_group import TileDBGroup

from typing import Optional, List, Dict


class TileDBArray(TileDBObject):
    """
    Wraps arrays from TileDB-Py by retaining a URI, verbose flag, etc.
    Also serves as an abstraction layer to hide TileDB-specific details from the API, unless
    requested.
    """

    def __init__(
        self,
        uri: str,
        name: str,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the TileDBObject constructor.
        """
        super().__init__(uri=uri, name=name, parent=parent)

    def _object_type(self):
        """
        This should be implemented by child classes and should return what tiledb.object_type(uri)
        returns for objects of a given type -- nominally 'group' or 'array'.
        """
        return "array"

    def _open(self):
        """
        Returns the TileDB array. The caller should do A.close() on the return value after finishing
        with it.
        """
        A = tiledb.open(self.uri)
        return A

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
        with tiledb.open(self.uri) as A:
            return A.schema

    def dim_names(self) -> List[str]:
        """
        Reads the dimension names from the schema: for example, ['obs_id', 'var_id'].
        """
        with tiledb.open(self.uri) as A:
            return [A.schema.domain.dim(i).name for i in range(A.schema.domain.ndim)]

    def dim_names_to_types(self) -> Dict[str, str]:
        """
        Returns a dict mapping from dimension name to dimension type.
        """
        with tiledb.open(self.uri) as A:
            dom = A.schema.domain
            return {dom.dim(i).name: dom.dim(i).dtype for i in range(dom.ndim)}

    def attr_names(self) -> List[str]:
        """
        Reads the attribute names from the schema: for example, the list of column names in a dataframe.
        """
        with tiledb.open(self.uri) as A:
            return [A.schema.attr(i).name for i in range(A.schema.nattr)]

    def attr_names_to_types(self) -> Dict[str, str]:
        """
        Returns a dict mapping from attribute name to attribute type.
        """
        with tiledb.open(self.uri) as A:
            schema = A.schema
            return {
                schema.attr(i).name: schema.attr(i).dtype for i in range(schema.nattr)
            }

    def has_attr_name(self, attr_name: str) -> bool:
        """
        Returns true if the array has the specified attribute name, false otherwise.
        """
        return attr_name in self.attr_names()
