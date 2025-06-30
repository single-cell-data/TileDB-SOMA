#!/usr/bin/env python3

import pytest

import tiledbsoma

from .common import TestReadPythonWriteR


class TestCharacterMetadataWriteRReadPython(TestReadPythonWriteR):
    @pytest.fixture(scope="class")
    def R_character(self):
        base_script = f"""
        library(tiledbsoma)

        exp <- SOMAExperimentCreate("{self.uri}")
        exp$close()
        """
        self.execute_R_script(base_script)

    def test_py_character(self, R_character):
        with tiledbsoma.open(self.uri) as exp:
            for key in exp.metadata:
                assert isinstance(exp.metadata.get(key), str)
