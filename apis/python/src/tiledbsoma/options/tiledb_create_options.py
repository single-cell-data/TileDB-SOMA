from typing import (
    Any,
    Dict,
    Iterator,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Type,
    Union,
    cast,
)

import attrs
import tiledb
from attrs import field
from somacore import options

DEFAULT_TILE_ORDER = "row-major"
DEFAULT_CELL_ORDER = "row-major"
DEFAULT_DATAFRAME_DIM_ZSTD_LEVEL = 3
DEFAULT_SPARSE_ND_ARRAY_DIM_ZSTD_LEVEL = 3
DEFAULT_WRITE_X_CHUNKED = True
DEFAULT_GOAL_CHUNK_NNZ = 200_000_000
DEFAULT_OFFSET_FILTERS = (
    "DoubleDeltaFilter",
    "BitWidthReductionFilter",
    "ZstdFilter",
)
DEFAULT_VALIDITY_FILTERS = None
DEFAULT_TILE_EXTENT = 2048
# TODO: pending further work on
#  https://github.com/single-cell-data/TileDB-SOMA/issues/27
# DEFAULT_OBS_EXTENT = 256
# DEFAULT_VAR_EXTENT = 2048
# DEFAULT_X_CAPACITY = 100000
# DEFAULT_MAX_THREAD_POOL_WORKERS = 8

StrOrMap = Union[str, Mapping[str, Any]]


@attrs.define(frozen=True)
class TileDBCreateOptions(Mapping[str, Any]):
    """
    Provides strongly-typed access and default values for `platform_config` options stored under the "tiledb"->"create"
    Mapping keys.

    Intended for internal use only.
    """

    _config: Mapping[str, Any] = field(factory=dict)

    @classmethod
    def from_platform_config(
        cls, platform_config: Optional[options.PlatformConfig] = None
    ) -> "TileDBCreateOptions":
        return TileDBCreateOptions(
            (platform_config or {}).get("tiledb", {}).get("create", {})
        )

    def dataframe_dim_zstd_level(self) -> int:
        return self.get("dataframe_dim_zstd_level", DEFAULT_DATAFRAME_DIM_ZSTD_LEVEL)

    def sparse_nd_array_dim_zstd_level(self) -> int:
        return self.get(
            "sparse_nd_array_dim_zstd_level", DEFAULT_SPARSE_ND_ARRAY_DIM_ZSTD_LEVEL
        )

    def write_X_chunked(self) -> bool:
        return self.get("write_X_chunked", DEFAULT_WRITE_X_CHUNKED)

    def goal_chunk_nnz(self) -> int:
        return self.get("goal_chunk_nnz", DEFAULT_GOAL_CHUNK_NNZ)

    def offsets_filters(
        self,
        default: Sequence[StrOrMap] = DEFAULT_OFFSET_FILTERS,
    ) -> Sequence[tiledb.Filter]:
        return _build_filters(self.get("offsets_filters", default))

    def validity_filters(
        self,
        default: Optional[Sequence[StrOrMap]] = DEFAULT_VALIDITY_FILTERS,
    ) -> Optional[Sequence[tiledb.Filter]]:
        items = self.get("validity_filters", default)
        if items is None:
            return None
        return _build_filters(items)

    def cell_tile_orders(self) -> Tuple[Optional[str], Optional[str]]:
        """Returns the cell and tile orders that should be used.

        If *neither* ``cell_order`` nor ``tile_order`` is present, only in this
        case will we use the default values provided.
        """
        if "cell_order" in self or "tile_order" in self:
            return self.get("cell_order"), self.get("tile_order")

        return DEFAULT_CELL_ORDER, DEFAULT_TILE_ORDER

    def dim_filters(
        self, dim: str, default: Sequence[StrOrMap] = ()
    ) -> Sequence[tiledb.Filter]:
        return _build_filters(self._dim(dim).get("filters", default))

    def dim_tile(self, dim: str, default: int = DEFAULT_TILE_EXTENT) -> int:
        return self._dim(dim).get("tile", default)

    def attr_filters(
        self, attr: str, default: Sequence[StrOrMap] = ()
    ) -> Sequence[tiledb.Filter]:
        return _build_filters(self._attr(attr).get("filters", default))

    def _dim(self, dim: str) -> Dict[str, Any]:
        return cast(Dict[str, Any], self.get("dims", {}).get(dim, {}))

    def _attr(self, attr: str) -> Dict[str, Any]:
        return cast(Dict[str, Any], self.get("attrs", {}).get(attr, {}))

    # Mapping implementation:

    def __getitem__(self, it: str) -> Any:
        return self._config[it]

    def __len__(self) -> int:
        return len(self._config)

    def __iter__(self) -> Iterator[str]:
        return iter(self._config)


_FILTERS: Mapping[str, Type[tiledb.Filter]] = {
    cls.__name__: cls for cls in tiledb.FilterList.filter_type_cc_to_python.values()
}


def _build_filters(items: Any) -> Sequence[tiledb.Filter]:
    return tuple(map(_build_filter, items))


def _build_filter(item: Any) -> tiledb.Filter:
    kwargs: Dict[str, Any]
    if isinstance(item, str):
        cls_name = item
        kwargs = {}
    elif isinstance(item, Mapping):
        kwargs = dict(item)
        try:
            cls_name = kwargs.pop("_type")
        except KeyError as ke:
            raise ValueError(
                "filter specification must include a `_type` key with the filter type"
            ) from ke
    else:
        raise TypeError(
            f"filter specification must be a string name or a dict, not {type(item)}"
        )

    try:
        cls = _FILTERS[cls_name]
    except KeyError as ke:
        raise ValueError(f"no filter named {cls_name!r}") from ke
    return cls(**kwargs)
