## Ingesting into SOMAs

A **SOMA** (Stack of Matrices, Annotated) is a [unified single-cell data model and API](https://github.com/single-cell-data/SOMA). A SOMA contains the kinds of data that belong in a single-cell dataset: `obs` and `var` matrices, `X` layers, and so on, offering **write-once, read-many semantics** via Python and R toolchains, with I/O to AnnData and Seurat formats, and interoperability with Scanpy and Seurat toolchains.

Use the [ingestor script](../tools/ingestor) to ingest [files you've downloaded](obtaining-data-files.md) into SOMAs:

```
tools/ingestor -o /mini-corpus/tiledb-data -n /mini-corpus/anndata/0cfab2d4-1b79-444e-8cbe-2ca9671ca85e.h5ad
tools/ingestor -o /mini-corpus/tiledb-data -n /mini-corpus/anndata/10x_pbmc68k_reduced.h5ad
...
```

Note this can take several hours total for multi-gigabyte datasets. The benefit of using an
optimized storage solution (with admittedly non-negligible ingest time) is that all subsequent
queries benefit from that optimized storage. In particular, various cross-corpus data queries shown
in these examples take just seconds or minutes.

A key point is **write once, read from multiple tools** -- in particular, using `tiledbsc-py` (this
package) or [`tiledbsc-r`](https://github.com/TileDB-Inc/tiledbsc) you can read SOMAs in either
language, regardless of which language was used to store them. This lets you use
best-in-class/state-of-the-art analysis algorithms, whichever language they're implemented in.

## Populate a SOMA collection

Use the [populator script](../tools/populate-soco) to mark these as members of a SOMA collection:

```
populate-soco -o /mini-corpus/soco -a /mini-corpus/tiledb-data/0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
populate-soco -o /mini-corpus/soco -a /mini-corpus/tiledb-data/10x_pbmc68k_reduced
...
```

or simply

```
populate-soco -o /mini-corpus/soco -a /mini-corpus/tiledb-data/*
```

Note this is quite quick.

As a keystroke-saver, use the `tools/ingestor` script's `--soco` option which will populate the SOMA
collection at ingest time, so you don't even have to run `populate-soco` as an afterstep.

```
tools/ingestor -o /mini-corpus/tiledb-data --soco -n /mini-corpus/anndata/0cfab2d4-1b79-444e-8cbe-2ca9671ca85e.h5ad
tools/ingestor -o /mini-corpus/tiledb-data --soco -n /mini-corpus/anndata/10x_pbmc68k_reduced.h5ad
```

## Names and URIs

Next let's start taking a look across the collection.

```
for soma in soco:
    print("%-40s %s" % (soma.name, soma.uri))
```

```
NAMES AND URIS
tabula-sapiens-immune                    file:///mini-corpus/tiledb-data/tabula-sapiens-immune
tabula-sapiens-epithelial                file:///mini-corpus/tiledb-data/tabula-sapiens-epithelial
integrated-human-lung-cell-atlas         file:///mini-corpus/tiledb-data/integrated-human-lung-cell-atlas
af9d8c03-696c-4997-bde8-8ef00844881b     file:///mini-corpus/tiledb-data/af9d8c03-696c-4997-bde8-8ef00844881b
subset_100_100                           file:///mini-corpus/tiledb-data/subset_100_100
pbmc3k_processed                         file:///mini-corpus/tiledb-data/pbmc3k_processed
pbmc3k-krilow                            file:///mini-corpus/tiledb-data/pbmc3k-krilow
issue-74                                 file:///mini-corpus/tiledb-data/issue-74
local2                                   file:///mini-corpus/tiledb-data/local2
developmental-single-cell-atlas-of-the-murine-lung file:///mini-corpus/tiledb-data/developmental-single-cell-atlas-of-the-murine-lung
single-cell-transcriptomes               file:///mini-corpus/tiledb-data/single-cell-transcriptomes
tabula-sapiens-stromal                   file:///mini-corpus/tiledb-data/tabula-sapiens-stromal
azimuth-meta-analysis                    file:///mini-corpus/tiledb-data/azimuth-meta-analysis
vieira19_Alveoli_and_parenchyma_anonymised.processed file:///mini-corpus/tiledb-data/vieira19_Alveoli_and_parenchyma_anonymised.processed
local3                                   file:///mini-corpus/tiledb-data/local3
Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection file:///mini-corpus/tiledb-data/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection
4056cbab-2a32-4c9e-a55f-c930bc793fb6     file:///mini-corpus/tiledb-data/4056cbab-2a32-4c9e-a55f-c930bc793fb6
longitudinal-profiling-49                file:///mini-corpus/tiledb-data/longitudinal-profiling-49
human-kidney-tumors-wilms                file:///mini-corpus/tiledb-data/human-kidney-tumors-wilms
issue-69                                 file:///mini-corpus/tiledb-data/issue-69
autoimmunity-pbmcs                       file:///mini-corpus/tiledb-data/autoimmunity-pbmcs
brown-adipose-tissue-mouse               file:///mini-corpus/tiledb-data/brown-adipose-tissue-mouse
d4db74ad-a129-4b1a-b9da-1b30db86bbe4-issue-74 file:///mini-corpus/tiledb-data/d4db74ad-a129-4b1a-b9da-1b30db86bbe4-issue-74
pbmc-small                               file:///mini-corpus/tiledb-data/pbmc-small
10x_pbmc68k_reduced                      file:///mini-corpus/tiledb-data/10x_pbmc68k_reduced
Puck_200903_10                           file:///mini-corpus/tiledb-data/Puck_200903_10
0cfab2d4-1b79-444e-8cbe-2ca9671ca85e     file:///mini-corpus/tiledb-data/0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
acute-covid19-cohort                     file:///mini-corpus/tiledb-data/acute-covid19-cohort
adult-mouse-cortical-cell-taxonomy       file:///mini-corpus/tiledb-data/adult-mouse-cortical-cell-taxonomy
```
