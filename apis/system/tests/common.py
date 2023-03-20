import subprocess


def to_R(x):
    if isinstance(x, bool):
        return str(x).upper()
    elif isinstance(x, str):
        return f'"{x}"'
    else:
        return str(x)


def create_R_list(xs):
    return ",".join([to_R(x) for x in xs])


class TestWritePythonReadR:
    @classmethod
    def setup_class(cls):
        import tempfile

        cls.td = tempfile.TemporaryDirectory()
        path = cls.td.name
        cls.uri = f"{path}/test.soma"

    @classmethod
    def teardown_class(cls):
        cls.td.cleanup()

    def r_assert(self, code: str):
        R_script = self.base_script() + code
        print(R_script)
        with open("test-dataframe-read.R", "w") as f:
            f.write(R_script)
        try:
            subprocess.run(["Rscript", "test-dataframe-read.R"], check=True)
        except subprocess.CalledProcessError:
            raise AssertionError(f"R assertion failed: {code}")
            # raise e
