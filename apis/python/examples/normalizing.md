## Normalizing a SOMA collection

The [cartographer script](../examples/cartographer.py) shows an example of how to take a collection
of H5AD files -- and/or already-ingested SOMAs -- and _normalize_ them into a uniform collection.

This is an alternative to using the [ingestor](../tools/ingestor) script -- the ingestor script pulls in data as-is, while this
cartographer does a strongly opinionated normalization. You might thinkof the cartographer script as a template for your own
organization-specific opinionated normalization.

```
examples/cartographer.py -v /Users/testuser/mini-corpus/atlas add-h5ad file-01.h5ad
examples/cartographer.py -v /Users/testuser/mini-corpus/atlas add-h5ad file-02.h5ad
...
```

```
examples/cartographer.py -v /Users/testuser/mini-corpus/atlas add-soma soma-01
examples/cartographer.py -v /Users/testuser/mini-corpus/atlas add-soma soma-02
...
```

After this, all the SOMAs in the collection will have the same schema:

```
>>> soco = tiledbsc.SOMACollection("/Users/testuser/mini-corpus/atlas")
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
