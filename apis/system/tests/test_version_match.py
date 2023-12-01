import tiledb

from .common import BasePythonRInterop


class TestVersionMatch(BasePythonRInterop):
    def test_version_match(self):
        """
        Verifies that the TileDB version of R and Python match. If they don't, the roundtrip
        testing will likely fail.
        """
        version = tiledb.libtiledb.version()
        major, minor = version[0], version[1]

        self.execute_R_script(
            f"""
        library("tiledb")
        version = tiledb_version()
        stopifnot(as.integer(version["major"]) == {major})
        stopifnot(as.integer(version["minor"]) >= {minor})
        """
        )
