from __future__ import annotations

import pathlib
import urllib.parse
from contextlib import contextmanager
from typing import Any

import tiledb.client


def parse_tiledb_uri(uri: str) -> tuple[str, str, str | None]:
    """Given a tiledb URI, return (workspace, teamspace, path)."""

    p_uri = urllib.parse.urlparse(uri)
    workspace = p_uri.netloc

    p_path = pathlib.Path(p_uri.path)
    if len(p_path.parts) < 2:
        raise ValueError(f"Not a tiledb URI - missing teamspace: {uri}")
    teamspace = p_path.parts[1]

    path = None if len(p_path.parts) == 2 else str(pathlib.PurePath(*p_path.parts[2:]))

    return workspace, teamspace, path


def join_tiledb_uri(workspace: str, teamspace: str, path: tuple[str] | str) -> str:
    path = (path,) if isinstance(path, str) else path
    return f"tiledb://{workspace}/{teamspace}/{'/'.join(path)}"


def get_asset_info(uri: str) -> dict[str, Any]:
    import tiledb.client

    _, teamspace, path = parse_tiledb_uri(uri)
    return tiledb.client.assets.get_asset(path, teamspace=teamspace).to_dict()


@contextmanager
def carrara_cleanup_asset(url: str) -> str:
    _, teamspace, path = parse_tiledb_uri(url)
    try:
        yield url
    finally:
        tiledb.client.assets.delete_asset(path, teamspace=teamspace, delete_storage=True)
