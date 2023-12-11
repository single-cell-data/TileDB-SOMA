# Issue:

# * tiledbsoma 1.5.1 with core 2.17 exists in GitHub, PyPI, r-universe, and Conda
# * tiledbsoma 1.6.0 with core 2.18 exists in GitHub, PyPI, r-universe, and Conda
# * We have a bugfix for a Python issue on release-1.6 and release-1.5
#   * main: https://github.com/single-cell-data/TileDB-SOMA/pull/1963
#   * release-1.6: https://github.com/single-cell-data/TileDB-SOMA/pull/1964 and this can go in a (not yet tagged) 1.6.1
#   * release-1.5: https://github.com/single-cell-data/TileDB-SOMA/pull/1968 and this can go in a (not yet tagged) 1.6.18
# * Problem:
#   * tiledb-r 0.22 (with core 2.18) now does exist
#   * The release-1.5 branch of TileDB-SOMA needs core 2.17
#   * The R DESCRIPTION language does not (syntactically) let us pin tiledb-r to == 0.21, only to >= 0.21
#   * CI fails with
#     * as.integer(version["minor"]) == 17 is not TRUE
# * Since 1.5.2 is a Python-only bugfix, we silence the CI alert here.

# import tiledb
#
# from .common import BasePythonRInterop
#
#
# class TestVersionMatch(BasePythonRInterop):
#    def test_version_match(self):
#        """
#        Verifies that the TileDB version of R and Python match. If they don't, the roundtrip
#        testing will likely fail.
#        """
#        version = tiledb.libtiledb.version()
#        major, minor = version[0], version[1]
#
#        self.execute_R_script(
#            f"""
#        library("tiledb")
#        version = tiledb_version()
#        stopifnot(as.integer(version["major"]) == {major})
#        stopifnot(as.integer(version["minor"]) == {minor})
#        """
#        )
