import json
import os
import subprocess
import sys
from typing import Dict

import psutil


def main():
    mem_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    physical_mem_gib = mem_bytes / (1024.0**3)

    def get_git_revision_hash() -> str:
        return (
            subprocess.check_output(["git", "rev-parse", "HEAD"])
            .decode("ascii")
            .strip()
        )

    stats: Dict[str, str] = {
        "uname": os.uname(),
        "total_virtual_mem": psutil.virtual_memory().total,
        "total_physical_mem": physical_mem_gib,
        "swap_mem": psutil.swap_memory(),
        "cpu_count": os.cpu_count(),
        "python_version": sys.version,
        "get_git_revision_hash": get_git_revision_hash(),
    }

    # This process is supposed to generate its output (context) in JSON format
    print(json.dumps(stats, indent=4))


if __name__ == "__main__":
    main()
