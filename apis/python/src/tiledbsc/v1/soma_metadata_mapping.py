from typing import Any, Dict, Iterator

# importing tiledbsc.v1.TileDBObject leads to a circular reference as TileDBObject imports us. This
# is, in turn, because this class requires a back-link to the underlying object -- hence,
# bidirectional references.
import tiledbsc.v1


class SOMAMetadataMapping:
    _underlying: "tiledbsc.v1.TileDBObject"

    def __init__(self, underlying: "tiledbsc.v1.TileDBObject"):
        self._underlying = underlying

    def get(self, key: str) -> Any:
        """
        Get the value associated with the key.
        """
        return self._get_all()[key]

    def has(self, key: str) -> bool:
        """
        Test for key existence.
        """
        return key in self._get_all()

    def set(self, key: str, value: Any) -> None:
        """
        Set the value associated with the key.
        """
        with self._underlying._tiledb_open("w") as M:
            M.meta[key] = value

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
        for k, v in self._get_all().items():
            # yield {k: v}
            yield (k, v)

    def __len__(self) -> int:
        """
        Get the length of the map, the number of keys present.
        """
        return len(self._get_all())

    def _get_all(self) -> Dict[str, Any]:
        """
        Retrieves the full metadata set from storage.
        """
        with self._underlying._tiledb_open("r") as M:
            return dict(M.meta)
