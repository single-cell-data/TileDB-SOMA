import tiledb
from .soma_options  import SOMAOptions
from .tiledb_object import TileDBObject
from .tiledb_group  import TileDBGroup

from typing import Optional

class TileDBArray(TileDBObject):
    """
    Wraps arrays from TileDB-Py by retaining a URI, verbose flag, etc.
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

    def object_type(self):
        """
        This should be implemented by child classes and should return what tiledb.object_type(uri)
        returns for objects of a given type -- nominally 'group' or 'array'.
        """
        return 'array'
