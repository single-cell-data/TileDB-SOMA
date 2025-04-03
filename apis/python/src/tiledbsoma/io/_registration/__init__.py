# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""
Support for soma_joinid remapping for append-mode ingestion.

This is an internal-use class; none of it is user-facing API. The user-facing API is
``tiledbsoma.io.register``.

The SOMA experiment ``obs``, ``var``, ``X``, etc. are indexed by soma_joinid.  Input AnnData/H5AD
``obs`` and ``var`` are indexed by an index column; ``X` et al. are indexed zero-up by row numbers
aligning with ``obs`` and ``var``. This raises the issue: if we have a SOMA experiment created
by ingesting two or more AnnData/H5AD files, how do we compute soma_joinid values for the full experiment?

Essential ideas:

- The input ``obs`` must have some (user-specified) column which contains values unique across all
  AnnData/H5AD inputs. Nominally this will be a cell barcode, and nominally this will be different
  for all input files. If one input has data for 100 cells and another has data for 200 cells,
  we will expect total input to have data for 300 cells. In particular, as multiple inputs
  are appended to a SOMA experiment, ``obs`` will grow taller.

- The input ``var`` must also have some (user-specified) column containing string identifiers.
  Nominally these are Ensembl IDs like ENSG00000142208, or HGNC IDs like AKT1. Nominally
  these will be the same for all input files, although it's to be expected that one input file
  may have data for some infrequently expressed genes that don't appear in other input files.
  In particular, as multiple inputs are appended to a SOMA experiment, ``var`` may gain
  a few rows here and there.

- Putting the last two together means that ``X`` (which is sparse) will mainly gain new
  rows (cells) for each input file, perhaps with some columns (genes) that didn't appear for
  previous input files.

The purpose of the registration mappings is to simply track a mapping from user-specified
``obs`` and ``var`` ID-column values to soma_joinid values.

Append-mode ingestion has a _registration pass_, which must be sequential and must use all input
files to compute join-ID mappings, and then an _ingestion pass_, which can be parallelized across
input files.

There are two kinds of mappings: _ambient mappings_ which contain string-to-int join-ID mappings
for _all_ inputs, and ID mappings which contain int-to-int offset-to-join-ID mappings for each _single_ input.

Example:

- Input 1 has obs IDs ``["AAAT", "ACTG", "AGAG"]`` numbered 0, 1, 2 within input 1.
- Input 1 has var IDs ``["AKT1", "APOE", "ESR1", "TP53", "VEGFA"]`` numbered 0, 1, 2, 3, 4 within input 1.
- Input 2 has obs IDs ``["CAAT", "CCTG", "CGAG"]`` numbered 0, 1, 2 within input 2.
- Input 2 has var IDs ``["APOE", "EGFR", "TP53", "VEGFA"]`` numbered 0, 1, 2, 3 within input 2.

Then registration produces ambient mappings like this:

- For ``obs``, ``AxisAmbientLabelMapping`` of

    AAAT:0
    ACTG:1
    AGAG:2
    CAAT:3
    CCTG:4
    CGAG:5

- For ``var``, ``AxisAmbientLabelMapping`` of

    AKT1:0
    APOE:1
    ESR1:2
    TP53:3
    VEGFA:4
    EGFR:5

- ``ExperimentAmbientLabelMapping`` containing these two axes.

This registration data is passed to the ingestor for each input file. Within the ingestion logic
itself and without user intervention, the ingestor selects out the int-to-int mappings from AnnData 0-up
offsets to registered SOMA join IDs like this:

- For input 1's ``obs``, ``AxisIDMapping`` of

    0:0 (for AAAT)
    1:1 (for ACTG)
    2:2 (for AGAG)

- For input 1's ``var``, ``AxisIDMapping`` of

    0:0 (for AKT1)
    1:1 (for APOE)
    2:2 (for ESR1)
    3:3 (for TP53)
    4:4 (for VEGFA)

- For input 2's ``obs``, ``AxisIDMapping`` of

    0:3 (for CAAT)
    1:4 (for CCTG)
    2:5 (for CGAG)

- For input 2's ``var``, ``AxisIDMapping`` of

    0:1 (for APOE)
    1:5 (for EGFR)
    2:3 (for TP53)
    3:4 (for VEGFA)
"""

from .ambient_label_mappings import (
    AxisAmbientLabelMapping,
    ExperimentAmbientLabelMapping,
)
from .id_mappings import AxisIDMapping, ExperimentIDMapping, get_dataframe_values

__all__ = (
    "AxisIDMapping",
    "AxisAmbientLabelMapping",
    "ExperimentIDMapping",
    "ExperimentAmbientLabelMapping",
    "get_dataframe_values",
)
