import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_amrfinderplus.feature_table import create_feature_table
from q2_amrfinderplus.types import AMRFinderPlusAnnotationsDirFmt


class TestFetchAMRFinderPlusDB(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    def test_create_feature_table(self):
        exp = pd.DataFrame(
            {
                "arsR": [0, 0, 0, 0, 1],
                "blaOXA": [0, 0, 1, 0, 0],
                "blaPDC": [0, 1, 0, 0, 0],
                "blaTEM": [0, 0, 0, 1, 0],
                "blaTEM-156": [1, 0, 0, 0, 0],
                "emrD3": [0, 0, 0, 0, 1],
            },
            index=["contig01", "contig02", "contig03", "contig08", "contig13"],
        )
        exp.index.name = "Contig id"
        exp.columns.name = "Gene symbol"
        annotations = AMRFinderPlusAnnotationsDirFmt(
            self.get_data_path("annotations_contigs"), mode="r"
        )
        obs = create_feature_table(annotations)
        pd.testing.assert_frame_equal(exp, obs)

    def test_value_error(self):
        annotations = AMRFinderPlusAnnotationsDirFmt(
            self.get_data_path("annotations_protein"), mode="r"
        )
        with self.assertRaisesRegex(ValueError, "solely from protein data"):
            create_feature_table(annotations)

    def test_empty_data_error(self):
        annotations = AMRFinderPlusAnnotationsDirFmt()
        with open(annotations.path / "sample1_amr_all_mutations.tsv", "w"):
            pass
        with self.assertRaisesRegex(ValueError, "File is empty"):
            create_feature_table(annotations)
