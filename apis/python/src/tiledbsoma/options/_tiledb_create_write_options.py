import warnings
from typing import (
    Any,
    Dict,
    Iterable,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Type,
    TypedDict,
    TypeVar,
    Union,
    cast,
)

import attrs as attrs_  # We use the name `attrs` later.
import attrs.validators as vld  # Short name because we use this a bunch.
from somacore import options
from typing_extensions import Self

import tiledb

from .._general_utilities import assert_version_before

# Most defaults are configured directly as default attribute values
# within TileDBCreateOptions.
DEFAULT_TILE_EXTENT = 2048
DEFAULT_CELL_ORDER = DEFAULT_TILE_ORDER = "row-major"
# TODO: pending further work on
#  https://github.com/single-cell-data/TileDB-SOMA/issues/27
# DEFAULT_X_CAPACITY = 100000
# DEFAULT_MAX_THREAD_POOL_WORKERS = 8

_DictFilterSpec = Mapping[str, object]
"""A format for specifying TileDB dimension/attribute filters and arguments.

The key ``_type`` is used as the name of the filter. Other entries in the
dictionary are passed as named arguments. For example,
``{"_type": "SomeFilter", "aggression": 5, "layers": 7}`` will call
``SomeFilter(aggression=5, layers=7)``."""
_FilterSpec = Union[str, _DictFilterSpec]
"""A declarative format for specifying filters:

- If a string, uses that filter with default settings. For instance,
  ``"OtherFilter"`` will call ``OtherFilter()``.
- If a dictionary, interprets as ``_DictFilterSpec``.
"""


class _DictColumnSpec(TypedDict, total=False):
    """Type specification for the dictionary used to configure a column."""

    filters: Union[Sequence[str], Mapping[str, Union[_FilterSpec]]]
    tile: int


# These functions have to appear first because they're used in the definition
# of TileDBCreateOptions.
def _normalize_filters(inputs: Iterable[_FilterSpec]) -> Tuple[_DictFilterSpec, ...]:
    if isinstance(inputs, str):
        raise TypeError(
            "filters must be a list of strings (or dicts), not a single string"
        )
    if not isinstance(inputs, Iterable):
        raise TypeError(
            f"filters must be a sequence of filter specs, not {type(inputs)}"
        )
    return tuple(_normalize_filter(spec) for spec in inputs)


# This exists because mypy does not currently (v1.3) support complex converters
# like converters.optional(inner_converter).
def _normalize_filters_optional(
    inputs: Optional[Iterable[_FilterSpec]],
) -> Optional[Tuple[_DictFilterSpec, ...]]:
    return None if inputs is None else _normalize_filters(inputs)


@attrs_.define(frozen=True, slots=True)
class _ColumnConfig:
    filters: Optional[Tuple[_DictFilterSpec, ...]] = attrs_.field(
        converter=_normalize_filters_optional
    )
    tile: Optional[int] = attrs_.field(validator=vld.optional(vld.instance_of(int)))

    @classmethod
    def from_dict(cls, input: _DictColumnSpec) -> Self:
        return cls(filters=input.get("filters"), tile=input.get("tile"))


def _normalize_columns(
    input: Mapping[str, _DictColumnSpec]
) -> Mapping[str, _ColumnConfig]:
    if not isinstance(input, Mapping):
        raise TypeError("column configuration must be a dictionary")
    return {
        col_name: _ColumnConfig.from_dict(value) for (col_name, value) in input.items()
    }


@attrs_.define(frozen=True, kw_only=True, slots=True)
class TileDBCreateOptions:
    """Tuning options used when creating new SOMA arrays.

    The attribute names of this object are identical to the keys expected in the
    cross-platform ``platform_config`` dict under the ``tiledb.create`` subkey,
    with values also accepted in that format (i.e., we create
    ``TileDBCreateOptions`` from a ``tiledb.create`` by directly calling
    ``TileDBCreateOptions(**tiledb_create_dict)``).
    """

    dataframe_dim_zstd_level: int = attrs_.field(
        validator=vld.instance_of(int), default=3
    )
    sparse_nd_array_dim_zstd_level: int = attrs_.field(
        validator=vld.instance_of(int), default=3
    )
    dense_nd_array_dim_zstd_level: int = attrs_.field(
        validator=vld.instance_of(int), default=3
    )
    write_X_chunked: bool = attrs_.field(validator=vld.instance_of(bool), default=True)
    goal_chunk_nnz: int = attrs_.field(
        validator=vld.instance_of(int), default=100_000_000
    )
    # We would prefer _remote_cap_nbytes as this is a server-side parameter
    # people should not be changing. However, leading underscores are not
    # accepted by the attrs framework.
    remote_cap_nbytes: int = attrs_.field(
        validator=vld.instance_of(int), default=2_400_000_000
    )
    capacity: int = attrs_.field(validator=vld.instance_of(int), default=100_000)
    offsets_filters: Tuple[_DictFilterSpec, ...] = attrs_.field(
        converter=_normalize_filters,
        default=(
            "DoubleDeltaFilter",
            "BitWidthReductionFilter",
            "ZstdFilter",
        ),
    )
    validity_filters: Optional[Tuple[_DictFilterSpec, ...]] = attrs_.field(
        converter=_normalize_filters_optional, default=None
    )
    allows_duplicates: bool = attrs_.field(
        validator=vld.instance_of(bool),
        default=False,
    )
    tile_order: Optional[str] = attrs_.field(
        validator=vld.optional(vld.instance_of(str)), default=None
    )
    cell_order: Optional[str] = attrs_.field(
        validator=vld.optional(vld.instance_of(str)), default=None
    )
    dims: Mapping[str, _ColumnConfig] = attrs_.field(
        factory=dict, converter=_normalize_columns
    )
    attrs: Mapping[str, _ColumnConfig] = attrs_.field(
        factory=dict, converter=_normalize_columns
    )

    @classmethod
    def from_platform_config(
        cls,
        platform_config: Union[
            options.PlatformConfig, "TileDBCreateOptions", None
        ] = None,
    ) -> Self:
        """Creates the object from a value passed in ``platform_config``.

        The value passed in should be the exact value passed into a public API
        method as the ``platform_config`` parameter. This function will extract
        the ``tiledb.create`` entry from a dict as needed.
        """
        create_entry = _dig_platform_config(platform_config, cls, ("tiledb", "create"))
        if isinstance(create_entry, dict):
            attrs: Tuple[attrs_.Attribute, ...] = cls.__attrs_attrs__  # type: ignore[type-arg]
            attr_names = frozenset(a.name for a in attrs)
            # Explicitly opt out of type-checking for these kwargs.
            filtered_create_entry: Dict[str, Any] = {
                key: value for (key, value) in create_entry.items() if key in attr_names
            }
            return cls(**filtered_create_entry)
        return create_entry

    def cell_tile_orders(self) -> Tuple[Optional[str], Optional[str]]:
        """Returns the cell and tile orders that should be used.

        If *neither* ``cell_order`` nor ``tile_order`` is present, only in this
        case will we use the default values provided.
        """
        if (self.cell_order, self.tile_order) == (None, None):
            return DEFAULT_CELL_ORDER, DEFAULT_TILE_ORDER
        return self.cell_order, self.tile_order

    def offsets_filters_tiledb(self) -> Tuple[tiledb.Filter, ...]:
        """Constructs the real TileDB Filters to use for offsets."""
        assert_version_before(1, 14)
        warnings.warn(
            "`offsets_filters_tiledb` is now deprecated for removal in 1.15 "
            "as we no longer support returning tiledb.Filter. "
            "Use `offsets_filters` instead.",
            DeprecationWarning,
        )

        return tuple(_build_filter(f) for f in self.offsets_filters)

    def validity_filters_tiledb(self) -> Optional[Tuple[tiledb.Filter, ...]]:
        """Constructs the real TileDB Filters to use for the validity map."""
        assert_version_before(1, 14)
        warnings.warn(
            "`validity_filters_tiledb` is now deprecated for removal in 1.15 "
            "as we no longer support returning tiledb.Filter. "
            "Use `validity_filters` instead.",
            DeprecationWarning,
        )
        if self.validity_filters is None:
            return None
        return tuple(_build_filter(f) for f in self.validity_filters)

    def dim_filters_tiledb(
        self, dim: str, default: Sequence[_FilterSpec] = ()
    ) -> Tuple[tiledb.Filter, ...]:
        """Constructs the real TileDB Filters to use for the named dimension."""
        assert_version_before(1, 14)
        warnings.warn(
            "`dim_filters_tiledb` is now deprecated for removal in 1.15 "
            "as we no longer support returning tiledb.Filter. "
            "Use `dims` instead.",
            DeprecationWarning,
        )
        return _filters_from(self.dims, dim, default)

    def dim_tile(self, dim_name: str, default: int = DEFAULT_TILE_EXTENT) -> int:
        """Returns the tile extent for the given dimension."""
        try:
            dim = self.dims[dim_name]
        except KeyError:
            return default
        return default if dim.tile is None else dim.tile

    def attr_filters_tiledb(
        self, name: str, default: Sequence[_FilterSpec] = ()
    ) -> Tuple[tiledb.Filter, ...]:
        """Constructs the real TileDB Filters to use for the named attribute."""
        assert_version_before(1, 14)
        warnings.warn(
            "`attr_filters_tiledb` is now deprecated for removal in 1.15 "
            "as we no longer support returning tiledb.Filter. "
            "Use `attrs` instead.",
            DeprecationWarning,
        )
        return _filters_from(self.attrs, name, default)


@attrs_.define(frozen=True, kw_only=True, slots=True)
class TileDBWriteOptions:
    """Tuning options used when writing to SOMA arrays."""

    sort_coords: bool = attrs_.field(validator=vld.instance_of(bool), default=True)
    consolidate_and_vacuum: Optional[bool] = attrs_.field(
        validator=vld.instance_of(bool), default=False
    )

    @classmethod
    def from_platform_config(
        cls,
        platform_config: Union[
            options.PlatformConfig, "TileDBWriteOptions", None
        ] = None,
    ) -> Self:
        """Creates the object from a value passed in ``platform_config``.

        The value passed in should be the exact value passed into a public API
        method as the ``platform_config`` parameter. This function will extract
        the ``tiledb.write`` entry from a dict as needed.
        """
        create_entry = _dig_platform_config(platform_config, cls, ("tiledb", "write"))
        if isinstance(create_entry, dict):
            attrs: Tuple[attrs_.Attribute, ...] = cls.__attrs_attrs__  # type: ignore[type-arg]
            attr_names = frozenset(a.name for a in attrs)
            # Explicitly opt out of type-checking for these kwargs.
            filered_create_entry: Dict[str, Any] = {
                key: value for (key, value) in create_entry.items() if key in attr_names
            }
            return cls(**filered_create_entry)
        return create_entry


_T = TypeVar("_T")


def _dig_platform_config(
    input: object, typ: Type[_T], full_path: Tuple[str, ...]
) -> Union[Dict[str, object], _T]:
    """Looks for an object of the given type in dictionaries.

    This is used to extract a valid object out of ``platform_config``. If an
    object of type ``typ`` is found, it is returned directly; otherwise we
    descend in the next dictionary key specified by ``full_path``. A ``dict``
    is allowed as the leaf element (to allow declarative configuration).
    """
    current = input
    path = full_path
    while path:  # Until we reach the leaf:
        if isinstance(current, typ):
            return current
        if not isinstance(current, dict):
            # Unrecognized type; return as an empty dict.
            return {}
        # Descend into the dict one more level.
        key, path = path[0], path[1:]
        try:
            current = current[key]
        except KeyError:
            # If the key isn't present, return as an empty dict.
            return {}
    # We've reached the leaf.  We validate this more closely.
    if not isinstance(current, (typ, dict)):
        path_dots = ".".join(full_path)
        raise TypeError(
            f"`{path_dots}` entry of `platform_config` must be"
            f" either a dict or `{typ.__name__}`, not {type(current)}"
        )
    # It's of the expected type! Return it.
    return current


#
# Filter handling and construction.
#

_FILTERS: Mapping[str, Type[tiledb.Filter]] = {
    cls.__name__: cls for cls in tiledb.FilterList.filter_type_cc_to_python.values()
}


def _normalize_filter(input: _FilterSpec) -> _DictFilterSpec:
    """Normalizes all filters to ``_DictFilterSpec`` format."""
    if isinstance(input, str):
        input = {"_type": input}
    if not isinstance(input, Mapping):
        raise TypeError(
            f"filters must be specified as a string or dict, not {type(input)}"
        )
    try:
        typ_name = input["_type"]
    except KeyError as ke:
        raise ValueError(
            "filter dicts must include a `_type` key with the filter name"
        ) from ke
    if not isinstance(typ_name, str):
        raise TypeError(f"filter name must be a str, not {type(typ_name)}")
    try:
        _FILTERS[typ_name]
    except KeyError as ke:
        raise ValueError(f"filter type {typ_name!r} unknown") from ke
    return dict(input)


def _filters_from(
    col_configs: Mapping[str, _ColumnConfig], name: str, default: Sequence[_FilterSpec]
) -> Tuple[tiledb.Filter, ...]:
    """Constructs the filters for the named column in ``col_configs``."""
    try:
        cfg = col_configs[name]
    except KeyError:
        maybe_filters = None
    else:
        maybe_filters = cfg.filters
    if maybe_filters is None:
        filters = _normalize_filters(default) or ()
    else:
        filters = maybe_filters
    return tuple(_build_filter(f) for f in filters)


def _build_filter(item: _DictFilterSpec) -> tiledb.Filter:
    """Build a single filter."""
    # Always make a copy here so we don't mutate the global state.
    # We have validated this earlier so we don't do extra checking here.
    kwargs = dict(item)
    cls_name = cast(str, kwargs.pop("_type"))
    cls = _FILTERS[cls_name]
    return cls(**kwargs)
