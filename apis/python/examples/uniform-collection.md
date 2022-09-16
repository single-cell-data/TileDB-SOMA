## Uniformizing a SOMA collection

The [uniformizer script](../examples/uniformizer.py) shows an example of how to take a collection
of H5AD files -- and/or already-ingested SOMAs -- and make them into a uniform collection.

:::{.callout-note}
This isn't necessary, or useful, for data exploration; you may find this information a helpful
guide if at any point in time your organization needs to construct an atlas.
:::

:::{.callout-note}
This is intended for bulk/batch/scripted jobs, as an alternative to using the
[ingestor](../tools/ingestor) script -- the ingestor script pulls in data as-is, while this
uniformizer is more strongly opinionated. You might think of this script as a template for your own
organization-specific opinionated uniformization.
:::

```
examples/uniformizer.py -v /Users/testuser/mini-corpus/atlas add-h5ad file-01.h5ad
examples/uniformizer.py -v /Users/testuser/mini-corpus/atlas add-h5ad file-02.h5ad
...
```

```
examples/uniformizer.py -v /Users/testuser/mini-corpus/atlas add-soma soma-01
examples/uniformizer.py -v /Users/testuser/mini-corpus/atlas add-soma soma-02
...
```

After this, all the SOMAs in the collection will have the same schema:

```
>>> soco = tiledbsoma.SOMACollection("/Users/testuser/mini-corpus/atlas")
>>> for soma in soco:
...     print(soma.obs.keys())
...
['assay_ontology_term_id', 'cell_type_ontology_term_id', 'development_stage_ontology_term_id', 'disease_ontology_term_id', 'ethnicity_ontology_term_id', 'organism_ontology_term_id', 'sex_ontology_term_id', 'tissue_ontology_term_id', 'assay', 'cell_type', 'development_stage', 'disease', 'ethnicity', 'organism', 'sex', 'tissue']
['assay_ontology_term_id', 'cell_type_ontology_term_id', 'development_stage_ontology_term_id', 'disease_ontology_term_id', 'ethnicity_ontology_term_id', 'organism_ontology_term_id', 'sex_ontology_term_id', 'tissue_ontology_term_id', 'assay', 'cell_type', 'development_stage', 'disease', 'ethnicity', 'organism', 'sex', 'tissue']
['assay_ontology_term_id', 'cell_type_ontology_term_id', 'development_stage_ontology_term_id', 'disease_ontology_term_id', 'ethnicity_ontology_term_id', 'organism_ontology_term_id', 'sex_ontology_term_id', 'tissue_ontology_term_id', 'assay', 'cell_type', 'development_stage', 'disease', 'ethnicity', 'organism', 'sex', 'tissue']
['assay_ontology_term_id', 'cell_type_ontology_term_id', 'development_stage_ontology_term_id', 'disease_ontology_term_id', 'ethnicity_ontology_term_id', 'organism_ontology_term_id', 'sex_ontology_term_id', 'tissue_ontology_term_id', 'assay', 'cell_type', 'development_stage', 'disease', 'ethnicity', 'organism', 'sex', 'tissue']
['assay_ontology_term_id', 'cell_type_ontology_term_id', 'development_stage_ontology_term_id', 'disease_ontology_term_id', 'ethnicity_ontology_term_id', 'organism_ontology_term_id', 'sex_ontology_term_id', 'tissue_ontology_term_id', 'assay', 'cell_type', 'development_stage', 'disease', 'ethnicity', 'organism', 'sex', 'tissue']
['assay_ontology_term_id', 'cell_type_ontology_term_id', 'development_stage_ontology_term_id', 'disease_ontology_term_id', 'ethnicity_ontology_term_id', 'organism_ontology_term_id', 'sex_ontology_term_id', 'tissue_ontology_term_id', 'assay', 'cell_type', 'development_stage', 'disease', 'ethnicity', 'organism', 'sex', 'tissue']
```
