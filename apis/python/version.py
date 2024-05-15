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
import subprocess
import sys
from datetime import date
from os.path import basename
from subprocess import CalledProcessError, check_output

RELEASE_VERSION_FILE = os.path.join(os.path.dirname(__file__), "RELEASE-VERSION")

# http://www.python.org/dev/peps/pep-0386/
_PEP386_SHORT_VERSION_RE = r"\d+(?:\.\d+)+(?:(?:[abc]|rc)\d+(?:\.\d+)*)?"
_PEP386_VERSION_RE = r"^%s(?:\.post\d+)?(?:\.dev\d+)?$" % _PEP386_SHORT_VERSION_RE
_GIT_DESCRIPTION_RE = r"^(?P<ver>%s)-(?P<commits>\d+)-g(?P<sha>[\da-f]+)$" % (
    _PEP386_SHORT_VERSION_RE
)


def err(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def readGitVersion():
    # NOTE: this will fail if on a fork with unsynchronized tags.
    #       use `git fetch --tags upstream`
    #       and `git push --tags <your fork>`
    try:
        proc = subprocess.Popen(
            ("git", "describe", "--long", "--tags", "--match", "[0-9]*.*"),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        data, stderr = proc.communicate()
        if proc.returncode:
            return None
        ver = data.decode().splitlines()[0].strip()
    except CalledProcessError:
        return None

    if not ver:
        return None
    m = re.search(_GIT_DESCRIPTION_RE, ver)
    if not m:
        err(
            "version: git description (%s) is invalid, " "ignoring\n" % ver,
        )
        return None

    commits = int(m.group("commits"))
    if not commits:
        return m.group("ver")
    else:
        return "%s.post%d.dev%d" % (m.group("ver"), commits, int(m.group("sha"), 16))


def readReleaseVersion():
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


def generateCalVersion():
    today = date.today().strftime("%Y.%m.%d")
    sha = check_output(["git", "log", "-1", "--format=%h"]).decode().rstrip("\n")
    sha_dec = int(sha, 16)
    return f"{today}.dev{sha_dec}"


def writeReleaseVersion(version):
    with open(RELEASE_VERSION_FILE, "w") as fd:
        print(version, file=fd)


def getVersion():
    release_version = readReleaseVersion()
    version = readGitVersion() or release_version
    if not version:
        version = generateCalVersion()
        err(
            f"No {basename(RELEASE_VERSION_FILE)} or Git version found, using calver {version}"
        )
    if version != release_version:
        writeReleaseVersion(version)
    return version
