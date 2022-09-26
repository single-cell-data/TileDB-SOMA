# Status

Temporary and experimental.

# Rationale

* R docs (currently at [https://github.com/TileDB-Inc/TileDB-SOMA](https://github.com/TileDB-Inc/TileDB-SOMA)) have a wonderful combination of API docs (generated from in-source-code doc-blocks) as well as hand-written long-form "vignette" material.
* For Python RST-style docs, I am not yet aware of a nice way to do that -- other than what's presented here.
* Tools like Sphinx and readthedocs are suitable for mapping a _single repo's single-language code-docs_ into a _single doc URL_. However, for this repo, we have Python API docs, Python examples/vignettes, and -- soon -- R docs as well. We wish to publish a _multi-lingual, multi-content doc tree_.

# Flow

* Source is in-source-code doc-blocks within `apis/python/src/tiledbsoma`, and hand-written long-form "vignette" material in `apis/python/examples`.
* The former are mapped to `.md` (intentionally not `.rst`) via [apis/python/mkmd.sh](apis/python/mkmd.sh). This requires `pydoc-markdown` already installed locally. (Nothing here in this initial experiment is CI-enabled at this point.)
* Then [Quarto](https://quarto.org) is used to map `.md` to `.html` via [_quarto.yml](_quarto.yml).
  * `quarto preview` for local preview.
  * `quarto render` to write static HTML into `docs/` which can then be published.
  * This `docs/` directory is artifacts-only and doesn't need to be committed to source control.
* Then this is synced to an AWS bucket which is used to serve static HTML content: [https://tiledb-singlecell-docs.s3.amazonaws.com/docs/overview.html](https://tiledb-singlecell-docs.s3.amazonaws.com/docs/overview.html).
  * [https://docs.aws.amazon.com/AmazonS3/latest/userguide/WebsiteAccessPermissionsReqd.html](https://docs.aws.amazon.com/AmazonS3/latest/userguide/WebsiteAccessPermissionsReqd.html)
