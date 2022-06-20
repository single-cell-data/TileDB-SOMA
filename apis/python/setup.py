import setuptools

if __name__ == "__main__":
    setuptools.setup(
        setup_requires=["setuptools_scm"],
        use_scm_version={
            "version_scheme": "guess-next-dev",
            "local_scheme": "dirty-tag",
            "write_to": "tiledb/ml/version.py",
        },
    )
