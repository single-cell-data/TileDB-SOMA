# TileDB-SOMA Development Policies

## Version policy

TileDB-SOMA uses a loose variation of [semantic versioning](https://semver.org/) to organize and govern API compatibility, deprecations and version numbers.

TileDB-SOMA also uses lifecycle tags to indicate the maturity of an interface. All public API should have a lifecycle tag in docstrings. These tags are patterned after the RStudio lifecycle stage model. Tags are:

- `experimental`: Under active development and may undergo significant and breaking changes.
- `maturing`: Under active development but the interface and behavior have stabilized and are unlikely to change significantly but breaking changes are still possible.
- `stable`: The interface is considered stable and breaking changes will be avoided where possible. Breaking changes that cannot be avoided will be accompanied by a major version bump.
- `deprecated`: The API is no longer recommended for use and will be removed in a future release.

If no tag is present, the state is `experimental`.

> [!NOTE]
> Prior to version 2.0, TileDB-SOMA had no explicit version policy. Versioning was most commonly `1.MAJOR.MINOR`,
> where breaking changes would occur in `1.X` releases, and any other change could occur on `1.X.Y` versions.
> Starting with `2.0.0`, versioning will be aligned with the SemVer `MAJOR.MINOR.PATCH` schema.

A release number is comprised of `MAJOR.MINOR.PATCH`.

For **non-experimental** interfaces and features:

- API breaking changes should occur only in `major` releases. These changes will be documented in the Python HISTORY.md or R NEWS.md change logs, with guidance on what has changed, and how to migrate code to the new version. When possible, advance warning will be provided prior to breaking changes.
- Deprecations may be introduced in `major` or `minor` releases. A deprecation will maintain existing functionality, but provide a warning about the upcoming API change, and any available guidance on migration. Deprecations must be documented in the appropriate change log.
- Warnings of future compatibility changes, which are not "deprecate and remove", may be introduced in a `major` or `minor` release. Future compatibility warnings must be included in the appropriate change log.
- A `patch` release may not introduce deprecations or future incompatibility change warnings.
- On disk format changes to the [TileDB storage format](https://github.com/TileDB-Inc/TileDB/tree/main/format_spec) will only occur on a `major` release.
- Changes to the [TileDB-SOMA encoding specification](https://github.com/single-cell-data/TileDB-SOMA/blob/main/encoding_specification.md) for non-experimental datatypes that are are either not backward-compatible or not forward-compatible (cannot be accessed from older versions of the TileDB-SOMA software) will only occur on a `major` release.


## Warning period

A best effort to use the following warning periods will be made. Any discrepancies will be highlighted in the changelog.

- Breaking API changes will provide notice of at least two `minor` releases or three months, whichever is greater, before effecting the change.
- Experimental features may be modified in a `major` or `minor` release without prior warning. However, where possible, advance warning of changes will be provided.
- Changes to the [TileDB-SOMA encoding specification](https://github.com/single-cell-data/TileDB-SOMA/blob/main/encoding_specification.md) that are not backward compatible will provide notice of at least two `minor` release or three months, whichever is greater, before effecting the change.
- No prior notice is given for changes to the [TileDB-SOMA encoding specification](https://github.com/single-cell-data/TileDB-SOMA/blob/main/encoding_specification.md) that are backward comaptible but not forward compatible.
- No prior notice is given for format changes to the [TileDB storage format](https://github.com/TileDB-Inc/TileDB/tree/main/format_spec).

## Implementing Deprecation

### Python deprecations

Upon deprecation of Python functionality:

- a `DeprecationWarning` must be issued by the deprecated functionality and a unit test provided to confirm warning issuance.
- affected lifecyle tags must be set to "deprecated", typically in the docstring and other user-visible documentation.
- docstrings must clearly indicate the scope of deprecation (e.g, function, a particular set of parameters, etc) and must include actionable guidance on migration.
- HISTORY.md must clearly indicate the deprecation and any guidance on migration.

For other breaking changes which are not "deprecate and remove":

- a `FutureWarning` must be issued by the deprecated functionality and a unit test provided to confirm warning issuance.
- docstrings must clearly indicate the scope of the breaking change (e.g, function, a particular set of parameters, etc) and must include actionable guidance on migration
- HISTORY.md must clearly indicate the breaking change and any guidance on migration.

Where entire classes, methods or functions are deprecated:

- utilize `typing_extensions.deprecated` decorator to provide both static and runtime warnings.
- update unit tests to ensure that the warning is issued (e.g., utilize `pytest.deprecated_call()`).
- ensure that `stacklevel` is set so that the warning comes from the invoking (caller's) code.

### R deprecations

To be defined.
