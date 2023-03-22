import subprocess


def embed_into_R_code(x):
    if isinstance(x, bool):
        return str(x).upper()
    elif isinstance(x, str):
        return f'"{x}"'
    else:
        return str(x)


def embed_python_list_into_R_code(xs):
    return ",".join([embed_into_R_code(x) for x in xs])


class BasePythonRInterop:
    @classmethod
    def setup_class(cls):
        import tempfile

        cls.tempdir = tempfile.TemporaryDirectory()
        path = cls.tempdir.name
        cls.uri = f"{path}/test.soma"

    @classmethod
    def teardown_class(cls):
        cls.tempdir.cleanup()

    def execute_R_script_file(self, script_name):
        subprocess.run(["Rscript", script_name], check=True)

    def execute_R_script(self, script):
        with open("test-dataframe-read.R", "w") as f:
            f.write(script)
        self.execute_R_script_file("test-dataframe-read.R")


class TestWritePythonReadR(BasePythonRInterop):
    def r_assert(self, code: str):
        R_script = self.base_R_script() + code
        try:
            self.execute_R_script(R_script)
        except subprocess.CalledProcessError:
            raise AssertionError(f"R assertion failed: {code}")
            # raise e


class TestReadPythonWriteR(BasePythonRInterop):
    pass
