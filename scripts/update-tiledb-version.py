#!/usr/bin/env python3

import argparse
import hashlib
import os
import re
from fileinput import FileInput
from subprocess import run
from urllib.request import urlopen


def hash_url_file(url):
    """Return SHA1 hash of the file located at the provided url."""

    BLOCK_SIZE = 65536
    hash = hashlib.sha1()
    with urlopen(url) as fp:
        while True:
            data = fp.read(BLOCK_SIZE)
            if not data:
                return hash.hexdigest()
            hash.update(data)


def get_version_hash(version):
    cmd = "git ls-remote --tags https://github.com/TileDB-Inc/TileDB.git"
    output = run(cmd, shell=True, capture_output=True).stdout.decode()

    m = re.search(rf"\s(\S+)\s+refs/tags/{version}\s", output)
    if m:
        return m.group(1)[0:7]

    print(output)
    print(f"Error: version {version} not found.")
    exit(1)


def main(args):
    old_version = None
    old_hash = None
    sha1 = None

    new_hash = get_version_hash(args.version)

    filepath = (
        f"{os.path.dirname(__file__)}/../libtiledbsc/cmake/Modules/FindTileDB_EP.cmake"
    )
    filepath = os.path.realpath(filepath)
    print(f"Updating {filepath}")
    print(f"  new version = {args.version}-{new_hash}")
    print(f"  computing SHA1 hashes...")

    # all "print" statements in this "with" block go to the "fp" file
    with FileInput(filepath, inplace=True) as fp:
        for line in fp:
            line = line.rstrip()

            if old_version is None:
                m = re.search(r"TileDB/releases/download/(.*?)/.*-(.*)\.zip", line)
                if m:
                    old_version = m.group(1)
                    old_hash = m.group(2)

            if old_version is not None:
                # modify url
                if "https://" in line:
                    line = line.replace(old_version, args.version)
                    line = line.replace(old_hash, new_hash)

                # update sha1 value computed on previous line
                if sha1 is not None:
                    if "URL_HASH" in line:
                        line = re.sub(r"SHA1=.*", f"SHA1={sha1}", line)
                    else:
                        line = re.sub(r'".*"', f'"{sha1}"', line)
                    sha1 = None
                else:
                    m = re.search(r'"(https://.*)"', line)
                    if m:
                        sha1 = hash_url_file(m.group(1))

            # print line to file
            print(line)

    print(f"  old version = {old_version}-{old_hash}")


if __name__ == "__main__":
    description = "Update FindTileDB_EP.cmake with a new TileDB version"
    epilog = f"Example: {__file__} 2.5.4"
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument("version", help="new TileDB version")
    args = parser.parse_args()

    main(args)
