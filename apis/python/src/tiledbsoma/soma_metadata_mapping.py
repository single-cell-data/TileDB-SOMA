from typing import Any, Dict, Iterator, MutableMapping, Union

# importing tiledbsoma.TileDBObject leads to a circular reference as TileDBObject imports us. This
# is, in turn, because this class requires a back-link to the underlying object -- hence,
# bidirectional references.
import tiledbsoma


class SOMAMetadataMapping(MutableMapping[str, Union[str, bool, int, float]]):
    _underlying: "tiledbsoma.TileDBObject"

    def __init__(self, underlying: "tiledbsoma.TileDBObject"):
        self._underlying = underlying

    def __delitem__(self, key: str) -> None:
        """
        Remove the key from the collection.
        """
        with self._underlying._tiledb_open("w") as M:
            del M.meta[key]

    def __iter__(self) -> Iterator[Any]:
        """
        Iterate over the collection.
        """
        return iter(self.as_dict())

    def __len__(self) -> int:
        """
        Get the length of the map, the number of keys present.
        """
        return len(self.as_dict())

    def as_dict(self) -> Dict[str, Any]:
        """
        Retrieves the full metadata set from storage.
        """
        with self._underlying._tiledb_open("r") as M:
            return dict(M.meta)

    def __getitem__(self, key: str) -> Any:
        """
        Implements ``item.metadata["key"]``.
        """
        with self._underlying._tiledb_open("r") as M:
            return M.meta[key]

    def __setitem__(self, key: str, value: Any) -> None:
        """
        Implements ``item.metadata["key"] = ...``.
        """
        with self._underlying._tiledb_open("w") as M:
            M.meta[key] = value
