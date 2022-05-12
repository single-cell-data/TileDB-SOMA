import tiledb
from .soma_options  import SOMAOptions
from .tiledb_object import TileDBObject

from typing import Optional, Union

class TileDBGroup(TileDBObject):
    """
    Wraps groups from TileDB-Py by retaining a URI, verbose flag, etc.
    """

    tiledb_group: Union[tiledb.Group, None]


    def __init__(
        self,
        uri: str,
        name: str,
        # Non-top-level objects can have a parent to propgate context, depth, etc.
        # What we really want to say is:
        # parent: Optional[TileDBGroup] = None,
        parent = None,
        # Top-level objects should specify these:
        soma_options: Optional[SOMAOptions] = None,
        verbose: Optional[bool] = True,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        See the TileDBObject constructor.
        """
        super().__init__(uri=uri, name=name, parent=parent)

        self.tiledb_group = None

    def object_type(self):
        """
        This should be implemented by child classes and should return what tiledb.object_type(uri)
        returns for objects of a given type -- nominally 'group' or 'array'.
        """
        return 'group'

    def create(self):
        """
        Creates the TileDB group data structure on disk/S3/cloud.
        """
        if self.tiledb_group != None:
            raise Exception("Attempt to create an already-open group")
        if self.verbose:
            print(f"{self.indent}Creating TileDB group {self.uri}")
        tiledb.group_create(uri=self.uri, ctx=self.ctx)

    def open(self, mode: str):
        """
        Mode must be "w" for write or "r" for read.
        """
        assert(mode in ['w', 'r'])
        if self.tiledb_group != None:
            raise Exception("Attempt to open an already-open group")
        if not self.exists():
            self.create()
        self.tiledb_group = tiledb.Group(self.uri, mode=mode, ctx=self.ctx)

    def close(self):
        """
        Should be done after open-with-write and add, or, open-with-read and read.
        """
        if self.tiledb_group == None:
            raise Exception("Attempt to close a non-open group")
        self.tiledb_group.close()
        self.tiledb_group = None

    def add(self, obj: TileDBObject):
        if self.tiledb_group == None:
            raise Exception("Attempt to write to a non-open group")
        if tiledb.object_type(obj.uri) is None:
            self.tiledb_group.add(uri=obj.uri, relative=False, name=obj.name)
        else:
            if self.verbose:
                print(f"{self.indent}Re-using existing group {obj.uri}")
