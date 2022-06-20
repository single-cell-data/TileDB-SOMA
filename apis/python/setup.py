import setuptools

if __name__ == "__main__":
    setuptools.setup(
        setup_requires=["setuptools_scm"],
        use_scm_version={
            "version_scheme": "guess-next-dev",
            "local_scheme": "dirty-tag",
            # This is weird. The write_to requires apis/python/..., as though the "pwd" is two
            # levels up from here.  Yet the root requires ../.., as though the "pwd" is right here.
            "write_to": "apis/python/version.py",
            "root": "../..",
        },
    )
