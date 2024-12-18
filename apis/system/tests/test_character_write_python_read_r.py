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
        return (
            """
        library(tiledbsoma)

        exp <- SOMAExperimentOpen("%s")
        md <- exp$get_metadata()
        for (i in seq_along(md)) {
          stopifnot(is.character(md[[i]]$name))
        }
        exp$close()
        """
            % self.uri
        )

    def test_r_character(self):
        self.r_assert("")
