## Normalizing a SOMA collection

The [cartographer script](../tools/cartographer) shows an example of how to take a collection
of H5AD files -- and/or already-ingested SOMAs -- and _normalize_ them into a uniform collection.

```
cartographer -v /Users/testuser/mini-corpus/atlas add-h5ad file-01.h5ad
cartographer -v /Users/testuser/mini-corpus/atlas add-h5ad file-02.h5ad
...
```

```
cartographer -v /Users/testuser/mini-corpus/atlas add-soma soma-01
cartographer -v /Users/testuser/mini-corpus/atlas add-soma soma-02
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
