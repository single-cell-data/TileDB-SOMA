from dataclasses import dataclass
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

import tiledb

from .types import PlatformConfig


@dataclass(frozen=True)
class TileDBPlatformConfig:
    """
    A place to put configuration options various users may wish to change.  These are mainly TileDB array-schema parameters.
    """

    # TODO: pending further work on
    # https://github.com/single-cell-data/TileDB-SOMA/issues/27
    obs_extent: int = 256
    var_extent: int = 2048
    X_capacity: int = 100000
    X_tile_order: str = "row-major"
    X_cell_order: str = "row-major"
    # https://github.com/single-cell-data/TileDB-SOMA/issues/27
    string_dim_zstd_level: int = 3
    write_X_chunked: bool = True
    goal_chunk_nnz: int = 200_000_000
    # Allows relocatability for local disk / S3, and correct behavior for TileDB Cloud
    member_uris_are_relative: Optional[bool] = None

    max_thread_pool_workers: int = 8


def from_param(pc: Optional[PlatformConfig]) -> "ParamWrapper":
    """Extracts the TileDB-specific platform config key."""
    pc = pc or {}
    return ParamWrapper(pc.get("tiledb", {}))


@dataclass(frozen=True)
class ParamWrapper(Mapping[str, Any]):
    """Wrapper for the TileDB key of the ``platform_config`` parameter.

    Provides some methods and proxies to the config dict itself.
    """

    _config: PlatformConfig

    def create_options(self) -> "CreateWrapper":
        return CreateWrapper(self.get("create", {}))

    # Mapping implementation:

    def __getitem__(self, it: str) -> Any:
        return self._config[it]

    def __len__(self) -> int:
        return len(self._config)

    def __iter__(self) -> Iterator[str]:
        return iter(self._config)


StrOrMap = Union[str, Mapping[str, Any]]


@dataclass(frozen=True)
class CreateWrapper(Mapping[str, Any]):

    _config: Mapping[str, Any]

    def offsets_filters(
        self,
        default: Sequence[StrOrMap] = (
            "DoubleDeltaFilter",
            "BitWidthReductionFilter",
            "ZstdFilter",
        ),
    ) -> Sequence[tiledb.Filter]:
        return _build_filters(self.get("offsets_filters", default))

    def cell_tile_orders(
        self,
        default_cell: Optional[str] = "row-major",
        default_tile: Optional[str] = "row-major",
    ) -> Tuple[Optional[str], Optional[str]]:
        """Returns the cell and tile orders that should be used.

        If *neither* ``cell_order`` nor ``tile_order`` is present, only in this
        case will we use the default values provided.
        """
        if "cell_order" in self or "tile_order" in self:
            return self.get("cell_order"), self.get("tile_order")

        return default_cell, default_tile

    def dim_filters(
        self, dim: str, default: Sequence[StrOrMap] = ()
    ) -> Sequence[tiledb.Filter]:
        return _build_filters(self._dim(dim).get("filters", default))

    def dim_tile(self, dim: str, default: int = 2048) -> int:
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
