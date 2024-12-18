#!/usr/bin/env python3

import pytest

import tiledbsoma

from .common import TestWritePythonReadR


class TestCharacterMetadataWritePythonReadR(TestWritePythonReadR):

    @pytest.fixture(scope="class")
    def experiment(self):
        exp = tiledbsoma.Experiment.create(self.uri)
        exp.close()

    def base_R_script(self):
        return f"""
        library(tiledbsoma)

        exp <- SOMAExperimentOpen("{self.uri}")
        md <- exp$get_metadata()
        """

    def test_r_character(self):
        self.r_assert(
            "stopifnot(all(vapply(md, \(x) is.character(x$name), logical(1L))))"
        )
