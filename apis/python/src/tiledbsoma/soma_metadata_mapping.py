from typing import Any, Dict, Iterator, Sequence

# importing tiledbsoma.TileDBObject leads to a circular reference as TileDBObject imports us. This
# is, in turn, because this class requires a back-link to the underlying object -- hence,
# bidirectional references.
import tiledbsoma


class SOMAMetadataMapping:
    _underlying: "tiledbsoma.TileDBObject"

    def __init__(self, underlying: "tiledbsoma.TileDBObject"):
        self._underlying = underlying

    def keys(self) -> Sequence[str]:
        """
        Returns the object's metadata keys as a list.
        """
        return list(self.items().keys())

    def get(self, key: str) -> Any:
        """
        Get the value associated with the key.
        """
        return self.items()[key]

    def has(self, key: str) -> bool:
        """
        Test for key existence.
        """
        return key in self.items()

    def set(self, key: str, value: Any) -> None:
        """
        Set the value associated with the key.
        """
        with self._underlying._tiledb_open("w") as M:
            M.meta[key] = value

    def __contains__(self, key: str) -> bool:
        """
        Answers whether the object's metadata contains the given key.
        """
        with self._underlying._tiledb_open("r") as M:
            return key in M.meta

    def __delete__(self, key: str) -> None:
        """
        Remove the key/value from the collection.
        """
        with self._underlying._tiledb_open("w") as M:
            del M.meta[key]

    def __iter__(self) -> Iterator[Any]:
        """
        Iterate over the collection.
        """
        for k, v in self.items().items():
            # yield {k: v}
            yield (k, v)

    def __len__(self) -> int:
        """
        Get the length of the map, the number of keys present.
        """
        return len(self.items())

    def items(self) -> Dict[str, Any]:
        """
        Retrieves the full metadata set from storage.
        """
        with self._underlying._tiledb_open("r") as M:
            return dict(M.meta)

    def __getitem__(self, key: str) -> Any:
        """
        Implements ``item.metadata["key"]``.
        """
        return self.get(key)

    def __setitem__(self, key: str, value: Any) -> None:
        """
        Implements ``item.metadata["key"] = ...``.
        """
        return self.set(key, value)
