import attrs
import h5py
import numpy as np
import pytest
import scipy.sparse as sp

ioutil = pytest.importorskip("tiledbsoma.io.spatial._util")


@attrs.define(frozen=True)
class SpaceRangerMatrixData:

    # Version
    version: tuple[int, int, int]

    # X
    nobs: int
    nvar: int
    data: np.ndarray[np.int32]
    obs_indices: np.ndarray[np.int64]
    var_indices: np.ndarray[np.int64]

    # obs (each should be length nobs)
    barcodes: np.ndarray[str]

    # var data (each should be length nvar)
    var_name: np.ndarray[str]
    gene_id: np.ndarray[str]
    feature_type: np.ndarray[str]
    genome: np.ndarray[str]

    # file
    filepath: str

    def __attrs_post_init__(self):
        with h5py.File(self.filepath, "w") as root:
            root.attrs["software_version"] = (
                f"spaceranger-{self.version[0]}.{self.version[1]}.{self.version[2]}"
            )

            # Create the matrix group.
            matrix_group = root.create_group("matrix")
            X = sp.csc_matrix(
                (self.data, (self.var_indices, self.obs_indices)),
                shape=np.array((self.nvar, self.nobs), np.int32),
            )
            matrix_group.create_dataset("barcodes", data=self.barcodes, dtype="S18")
            matrix_group.create_dataset("data", data=X.data, dtype=np.int32)
            matrix_group.create_dataset("indices", data=X.indices, dtype=np.int64)
            matrix_group.create_dataset("indptr", data=X.indptr, dtype=np.int64)
            matrix_group.create_dataset(
                "shape", data=np.array([self.nvar, self.nobs], dtype=np.int32)
            )

            # Create the feature group.
            features_group = matrix_group.create_group("features")
            features_group.create_dataset(
                "feature_type", data=self.feature_type, dtype="S15"
            )
            features_group.create_dataset("genome", data=self.genome, dtype="S6")
            features_group.create_dataset("id", data=self.gene_id, dtype="S15")
            features_group.create_dataset("name", data=self.var_name, dtype="S17")


@pytest.fixture(scope="session")
def fake_space_ranger_matrix_1(tmp_path_factory):
    # Create filepath.
    h5path = tmp_path_factory.mktemp("space_ranger") / "fake_space_ranger_matrix_1.h5"

    # Create test data.
    test_data = SpaceRangerMatrixData(
        version=(2, 0, 1),
        nobs=9,
        nvar=7,
        data=np.array([3, 8, 2, 4, 10, 3], dtype=np.int32),
        obs_indices=np.array([1, 2, 4, 4, 4, 6], dtype=np.int64),
        var_indices=np.array([0, 0, 3, 4, 5, 1], dtype=np.int64),
        barcodes=np.array(
            [
                "AAACTAACGTGGCGAC-1",
                "AAACTCGGTTCGCAAT-1",
                "AAACTCGTGATATAAG-1",
                "AAACTGCTGGCTCCAA-1",
                "AAACTTAATTGCACGC-1",
                "AAACTTGCAAACGTAT-1",
                "AAAGACATGAAGTTTA-1",
                "AAAGACCCAAGTCGCG-1",
                "AAAGACTGGGCGCTTT-1",
            ],
            dtype="S18",
        ),
        var_name=np.array(
            [
                "ISG15",
                "AL645608.1",
                "AGRN",
                "AL645608.5",
                "AL645608.8",
                "RNF223",
                "C1orf159",
            ],
            dtype="S15",
        ),
        gene_id=np.array(
            [
                "ENSG00000187608",
                "ENSG00000224969",
                "ENSG00000188157",
                "ENSG00000242590",
                "ENSG00000273443",
                "ENSG00000237330",
                "ENSG00000131591",
            ],
            dtype="S15",
        ),
        feature_type=np.array(7 * ["Gene Expression"], dtype="S15"),
        genome=np.array(
            [
                "ISG15",
                "AL645608.1",
                "AGRN",
                "AL645608.5",
                "AL645608.8",
                "RNF223",
                "C1orf159",
            ],
            dtype="S6",
        ),
        filepath=h5path,
    )
    return test_data


def check_reader(reader, data):
    """Check the data read by the :class:`TenXCountMatrixReader` against
    :class:`SpaceRangerMatrixData`.

    Either the data must already be loaded or :class:`SpaceRangerMatrixReader`
    must be open.
    """
    assert reader.software_version == data.version

    # Check shape.
    assert reader.nobs == data.nobs
    assert reader.nvar == data.nvar

    # Check X matrix.
    actual_X = sp.coo_matrix(
        (reader.data, (reader.obs_indices, reader.var_indices)),
        shape=(reader.nobs, reader.nvar),
    )
    expected_X = sp.coo_matrix(
        (data.data, (data.obs_indices, data.var_indices)),
        shape=(data.nobs, data.nvar),
    )
    np.testing.assert_equal(actual_X.todense(), expected_X.todense())

    # Check obs data.
    np.testing.assert_equal(reader.obs_id, data.barcodes.astype(str))

    # Check var data.
    np.testing.assert_equal(reader.var_id, data.var_name.astype(str))
    np.testing.assert_equal(reader.gene_id, data.gene_id.astype(str))
    np.testing.assert_equal(reader.feature_type, data.feature_type.astype(str))
    np.testing.assert_equal(reader.genome, data.genome.astype(str))

    # Check unique obs indices.
    actual_unique_obs_ind = np.sort(reader.unique_obs_indices().to_numpy())
    expected_unique_obs_ind = np.sort(np.unique(data.obs_indices))
    np.testing.assert_equal(actual_unique_obs_ind, expected_unique_obs_ind)

    # Check unique var indices.
    actual_unique_var_ind = np.sort(reader.unique_var_indices().to_numpy())
    expected_unique_var_ind = np.sort(np.unique(data.var_indices))
    np.testing.assert_equal(actual_unique_var_ind, expected_unique_var_ind)


def test_space_ranger_matrix_reader_direct_load(fake_space_ranger_matrix_1):

    # Open the reader, load the version and other data, and close.
    reader = ioutil.TenXCountMatrixReader(fake_space_ranger_matrix_1.filepath)
    reader.open()
    reader.software_version
    reader.load()
    reader.close()

    check_reader(reader, fake_space_ranger_matrix_1)


def test_space_ranger_matrix_reader_lazy_load(fake_space_ranger_matrix_1):

    with ioutil.TenXCountMatrixReader(fake_space_ranger_matrix_1.filepath) as reader:
        check_reader(reader, fake_space_ranger_matrix_1)
