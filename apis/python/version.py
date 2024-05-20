# -*- coding: utf-8 -*-

"""
NOTICE (mlin 2023-03-06): this script derives the Python package version number
based on the git history/tags. It also stores that info in a RELEASE-VERSION
file for use once source code has been distributed separately from the git repo
(e.g. in an sdist tarball). The script was originally obtained from:
  https://gist.github.com/mina86/8782771

It'd be preferable to use setuptools_scm instead of this ad hoc script.
Unfortunately, as of this writing setuptools_scm has issues with our repo
structure in which the Python package resides in a subdirectory as opposed to
the repo root. When setuptools builds the sdist tarball, it naturally homes
the Python package in the root. We haven't found a way to make setuptools_scm
handle the package location differing based on whether we're operating from the
source repo or an extracted sdist. The following setuptools_scm issues capture
this:
  https://github.com/pypa/setuptools_scm/issues/188
  https://github.com/pypa/setuptools_scm/issues/788

Note the problem doesn't necessarily become apparent until testing the sdist.

----

Calculates the current version number.

If possible, uses output of “git describe” modified to conform to the
versioning scheme that setuptools uses (see PEP 386).  Releases must be
labelled with annotated tags (signed tags are annotated) of the following
format:

   v<num>(.<num>)+ [ {a|b|c|rc} <num> (.<num>)* ]

If “git describe” returns an error (likely because we're in an unpacked copy
of a release tarball, rather than a git working copy), or returns a tag that
does not match the above format, version is read from RELEASE-VERSION file.

To use this script, simply import it your setup.py file, and use the results
of getVersion() as your package version:

    import version
    setup(
        version=version.getVersion(),
        .
        .
        .
    )

This will automatically update the RELEASE-VERSION file.  The RELEASE-VERSION
file should *not* be checked into git but it *should* be included in sdist
tarballs (as should version.py file).  To do this, run:

    echo include RELEASE-VERSION version.py >>MANIFEST.in
    echo RELEASE-VERSION >>.gitignore

With that setup, a new release can be labelled by simply invoking:

    git tag -s v1.0
"""

__author__ = (
    "Douglas Creager <dcreager@dcreager.net>",
    "Michal Nazarewicz <mina86@mina86.com>",
)
__license__ = "This file is placed into the public domain."
__maintainer__ = "Michal Nazarewicz"
__email__ = "mina86@mina86.com"

__all__ = "getVersion"


import os
import re
import shlex
import sys
from datetime import date
from os.path import basename
from subprocess import DEVNULL, CalledProcessError, check_output
from typing import List, Optional

RELEASE_VERSION_FILE = os.path.join(os.path.dirname(__file__), "RELEASE-VERSION")

# http://www.python.org/dev/peps/pep-0386/
_PEP386_SHORT_VERSION_RE = r"\d+(?:\.\d+)+(?:(?:[abc]|rc)\d+(?:\.\d+)*)?"
_PEP386_VERSION_RE = r"^%s(?:\.post\d+)?(?:\.dev\d+)?$" % _PEP386_SHORT_VERSION_RE
_GIT_DESCRIPTION_RE = (
    r"^(?P<ver>%s)-(?P<commits>\d+)-g(?P<sha>[\da-f]+)$" % _PEP386_SHORT_VERSION_RE
)


def err(*args, **kwargs):
    """Print to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def lines(*cmd, drop_trailing_newline: bool = True, **kwargs) -> List[str]:
    """Run a command and return its output as a list of lines.

    Strip trailing newlines, and drop the last line if it's empty, by default."""
    lns = [ln.rstrip("\n") for ln in check_output(cmd, **kwargs).decode().splitlines()]
    if lns and drop_trailing_newline and not lns[-1]:
        lns.pop()
    return lns


def line(*cmd, **kwargs) -> str:
    """Run a command, verify exactly one line of stdout, return it."""
    lns = lines(*cmd, **kwargs)
    if len(lns) != 1:
        raise RuntimeError(f"Expected 1 line, found {len(lns)}: {shlex.join(cmd)}")
    return lns[0]


def get_latest_tag() -> Optional[str]:
    """Return the most recent local Git tag of the form `[0-9].*.*` (or `None` if none exist)."""
    try:
        tags = lines("git", "tag", "--list", "--sort=v:refname", "[0-9].*.*")
        return tags[-1] if tags else None
    except CalledProcessError:
        return None


def get_latest_remote_tag(remote: str) -> str:
    """Return the most recent Git tag of the form `[0-9].*.*`, from a remote Git repository."""
    tags = lines("git", "ls-remote", "--tags", "--sort=v:refname", remote, "[0-9].*.*")
    if not tags:
        raise RuntimeError(f"No tags found in remote {remote}")
    return tags[-1].split(" ")[-1].split("/")[-1]


def get_sha_base10() -> int:
    """Return the current Git SHA, abbreviated and then converted to base 10.

    This is unfortunately necessary because PEP440 prohibits hexadecimal characters"""
    sha = line("git", "log", "-1", "--format=%h")
    return int(sha, 16)


def get_only_remote() -> Optional[str]:
    """Find the only remote Git repository, if one exists."""
    try:
        remotes = lines("git", "remote")
        if len(remotes) == 1:
            return remotes[0]
    except CalledProcessError:
        pass
    return None


def get_default_remote() -> Optional[str]:
    """Find a Git remote to parse a most recent release tag from.

    - If the current branch tracks a remote branch, use that remote
    - Otherwise, if there's only one remote, use that
    - Otherwise, return `None`
    """
    try:
        tracked_branch = line(
            "git", "rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"
        )
        tracked_remote = tracked_branch.split("/")[0]
        err(f"Parsed tracked remote {tracked_remote} from branch {tracked_branch}")
        return tracked_remote
    except CalledProcessError:
        remote = get_only_remote()
        if remote:
            err(f"Checking tags at default/only remote {remote}")
            return remote
        else:
            return None


def get_git_version() -> Optional[str]:
    """Construct a PEP440-compatible version string that encodes various Git state.

    - If `git describe` returns a plain release tag, use that.
    - Otherwise, it will return something like `1.10.2-5-gbabb931f2`, which we'd convert to
      `1.10.2.post5.dev50125681138` (abbreviated Git SHA gets converted to base 10, for
      PEP440-compliance).
    - However, if the `git describe` version starts with `1.5.0`, we do something else. 1.5.0 was
      the last release tag before we moved to release branches, so it ends up being returned for
      everything on the `main` branch. Instead:
        - Find the latest release tag in the local repo (or tracked remote, if there are no local
          tags).
        - Build a version string from it, e.g. `1.11.1.post0.dev61976836339` (again using the
          abbreviated Git SHA, converted to base 10 for PEP440 compliance).
    """
    try:
        git_version = line(
            "git", "describe", "--long", "--tags", "--match", "[0-9]*.*", stderr=DEVNULL
        )
    except CalledProcessError:
        git_version = None

    m = re.search(_GIT_DESCRIPTION_RE, git_version) if git_version else None
    ver = m.group("ver") if m else None

    # `1.5.0` (Nov '23) is typically the most recent tag that's an ancestor of `main`; subsequent
    # release tags all exist on release branches (by design).
    #
    # If `git describe` above returned `1.5.0` as the nearest tagged ancestor, synthesize a
    # more meaningful version number below:
    #
    # 1. Find the latest release tag in the local repo (or tracked remote, if there are no local
    #    tags, e.g. in case of a shallow clone).
    # 2. Return a PEP440-compatible version of the form `A.B.C.post0.devN`, where:
    #   - `A.B.C` is the most recent release tag in the repo, and
    #   - `N` is the current short Git SHA, converted to base 10.
    if not ver or ver.startswith("1.5.0"):
        latest_tag = get_latest_tag()
        if latest_tag:
            err(f"Git traversal returned {ver}, using latest local tag {latest_tag}")
        else:
            # Didn't find a suitable local tag, look for a tracked/default "remote", and find its
            # latest release tag
            tracked_remote = get_default_remote()
            if tracked_remote:
                try:
                    latest_tag = get_latest_remote_tag(tracked_remote)
                    err(
                        f"Git traversal returned {ver}, using latest tag {latest_tag} from tracked remote {tracked_remote}"
                    )
                except CalledProcessError:
                    err(f"Failed to find tags in remote {tracked_remote}")
                    return None
            else:
                err("Failed to find a suitable remote for tag traversal")
                return None

        try:
            sha_base10 = get_sha_base10()
            return f"{latest_tag}.post0.dev{sha_base10}"
        except CalledProcessError:
            err("Failed to find current SHA")
            return None
    else:
        commits = int(m.group("commits"))
        if commits:
            sha_base10 = int(m.group("sha"), 16)
            return f"{ver}.post{commits}.dev{sha_base10}"
        else:
            return ver


def read_release_version() -> Optional[str]:
    try:
        with open(RELEASE_VERSION_FILE) as fd:
            ver = fd.readline().strip()
        if not re.search(_PEP386_VERSION_RE, ver):
            err(
                "version: release version (%s) is invalid, "
                "will use it anyway\n" % ver,
            )
        return ver
    except FileNotFoundError:
        return None


def generate_cal_version():
    today = date.today().strftime("%Y.%m.%d")
    return f"{today}.dev{get_sha_base10()}"


def write_release_version(version):
    with open(RELEASE_VERSION_FILE, "w") as fd:
        print(version, file=fd)


def get_version():
    release_version = read_release_version()
    version = get_git_version()
    if not version:
        version = release_version
    if not version:
        version = generate_cal_version()
        err(
            f"No {basename(RELEASE_VERSION_FILE)} or Git version found, using calver {version}"
        )
    if version != release_version:
        write_release_version(version)
    return version


if __name__ == "__main__":
    print(get_version())
