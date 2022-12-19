from typing import Any, Dict, Iterator, MutableMapping

# importing tiledbsoma.TileDBObject leads to a circular reference as TileDBObject imports us. This
# is, in turn, because this class requires a back-link to the underlying object -- hence,
# bidirectional references.
import tiledbsoma


class MetadataMapping(MutableMapping[str, Any]):
    _underlying: "tiledbsoma.TileDBObject"

    def __init__(self, underlying: "tiledbsoma.TileDBObject"):
        self._underlying = underlying

    def __delitem__(self, key: str) -> None:
        """
        Remove the key from the metadata collection.
        """
        with self._underlying._tiledb_open("w") as M:
            del M.meta[key]

    def __iter__(self) -> Iterator[str]:
        """
        Return an iterator over the metadata collection.
        """
        return iter(self.as_dict())

    def __len__(self) -> int:
        """
        Return the number of elements in the metadata collection.

        Returns
        -------
        int
            The number of elements in the metadata collection.
        """
        return len(self.as_dict())

    def as_dict(self) -> Dict[str, Any]:
        """
        Return the metadata collection as a ``dict``.

        Returns
        -------
        dict[str, any]
            The contents of the metadata collection.
        """
        with self._underlying._tiledb_open("r") as M:
            return dict(M.meta)

    def __getitem__(self, key: str) -> Any:
        """
        Return the metadata element specified by ``key``.

        Parameters
        ----------
        key : str
            The name of the item.
        """
        with self._underlying._tiledb_open("r") as M:
            return M.meta[key]

    def __setitem__(self, key: str, value: Any) -> None:
        """
        Update the metadata collection with a new element.

        Parameters
        ----------
        key : str
            The metadata element name.
        value : any
            The metadata element value. Must be a primitive type (int,
            float, bool) or string.

        Returns
        -------
        None
        """

        if type(key) != str:
            raise TypeError("Metadata keys must be a string.")

        with self._underlying._tiledb_open("w") as M:
            M.meta[key] = value
