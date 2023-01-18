from typing import TypedDict

from tiledbsoma.config import TileDBCreateOptions


class TileDBPlatformConfig(TypedDict):
    create: TileDBCreateOptions


# complies with somacore.options.PlatformConfig
class PlatformConfig(TypedDict):
    tiledb: TileDBPlatformConfig
