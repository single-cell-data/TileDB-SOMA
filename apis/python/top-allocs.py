#!/usr/bin/env python

import json
import re
from importlib.metadata import version

from click import argument, command, option
from utz import iec, o, proc

REPOS = {
    "dask": "dask/dask",
    "cloudpickle": "cloudpipe/cloudpickle",
    "msgpack": "msgpack/msgpack-python",
    "tiledbsoma": "single-cell-data/TileDB-SOMA",
}


@command
@option("-m", "--markdown", is_flag=True)
@argument("stats_path")
def main(markdown: bool, stats_path: str):
    with open(stats_path) as f:
        stats = o(json.load(f))

    if markdown:
        print("<div>")
        print()
    for alloc in stats.top_allocs:
        loc = alloc["location"]
        for k in ["site-packages/", "python/src/"]:
            if k in loc:
                i = loc.index(k)
                loc = loc[i + len(k) :]

        if loc != "<stack trace unavailable>":
            loc, line = loc.rsplit(":", 1)
            line = int(line)
            pkg, _ = loc.split("/", 1)
            repo = REPOS[pkg]
            if pkg == "tiledbsoma":
                path = "apis/python/src/"
                v = proc.line("git", "log", "-1", "--format=%h", log=False)
            else:
                path = ""
                v = version(pkg)
                if pkg != "dask":
                    v = f"v{v}"
            txt = re.sub("_", r"\_", loc)
            url = f"http://github.com/{repo}/blob/{v}/{path}{loc}#L{line}"
            sz_str = f"{iec(alloc['size'])}"
            path_str = f"[{txt}]({url})"
        else:
            sz_str = f"{iec(alloc['size'])}"
            path_str = f"`{loc}`"
        if markdown:
            print(f"- {path_str}: {sz_str}")
        else:
            print(f"{sz_str}: {path_str}")
    if markdown:
        print("</div>")


if __name__ == "__main__":
    main()
