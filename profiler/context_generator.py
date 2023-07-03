import json
import os
import subprocess
import sys
from typing import Dict

import psutil


def host_context() -> Dict[str, str]:
    physical_mem_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")

    def get_git_revision_hash() -> str:
        return (
            subprocess.check_output(["git", "rev-parse", "HEAD"])
            .decode("ascii")
            .strip()
        )

    swap_mem = psutil.swap_memory()

    return {
        "uname": os.uname(),
        "total_virtual_mem_bytes": psutil.virtual_memory().total,
        "total_physical_mem_bytes": physical_mem_bytes,
        "swap_mem_total_bytes": swap_mem.total,
        "swap_mem_used_bytes": swap_mem.used,
        "swap_mem_free_bytes": swap_mem.free,
        "swap_mem_pct_bytes": swap_mem.percent,
        "swap_mem_sin_bytes": swap_mem.sin,
        "swap_mem_sout_bytes": swap_mem.sout,
        "cpu_count": os.cpu_count(),
        "python_version": sys.version_info,
        "get_git_revision_hash": get_git_revision_hash(),
    }


if __name__ == "__main__":
    print(json.dumps(host_context()))
