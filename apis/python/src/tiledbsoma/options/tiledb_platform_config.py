from typing import TypedDict

from tiledbsoma.options import TileDBCreateOptions


class TileDBPlatformConfig(TypedDict):
    create: TileDBCreateOptions
