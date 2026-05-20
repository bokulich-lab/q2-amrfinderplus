import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_amrfinderplus.feature_table import create_feature_table
from q2_amrfinderplus.types import AMRFinderPlusAnnotationsDirFmt


class TestFetchAMRFinderPlusDB(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    def setUp(self):
        super().setUp()
        self.exp_gene_table = pd.DataFrame(
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
        self.exp_gene_table.index.name = "Contig id"

    def test_create_feature_table_gene(self):
        exp = self.exp_gene_table.copy()
        exp.columns.name = "Gene symbol"
        annotations = AMRFinderPlusAnnotationsDirFmt(
            self.get_data_path("annotations_contigs_1"), mode="r"
        )
        obs = create_feature_table(annotations, level="gene")
        pd.testing.assert_frame_equal(exp, obs)

    def test_create_feature_table_gene_new_header(self):
        exp = self.exp_gene_table.copy()
        exp.columns.name = "Element symbol"
        annotations = AMRFinderPlusAnnotationsDirFmt(
            self.get_data_path("annotations_contigs_new_header"), mode="r"
        )
        obs = create_feature_table(annotations, level="gene")
        pd.testing.assert_frame_equal(exp, obs)

    def test_create_feature_table_class(self):
        exp = pd.DataFrame(
            {
                "ARSENIC": [0, 0, 0, 0, 1],
                "BETA-LACTAM": [1, 1, 1, 1, 0],
                "EFFLUX": [0, 0, 0, 0, 1],
            },
            index=["contig01", "contig02", "contig03", "contig08", "contig13"],
        )
        exp.index.name = "Contig id"
        exp.columns.name = "Class"
        annotations = AMRFinderPlusAnnotationsDirFmt(
            self.get_data_path("annotations_contigs_1"), mode="r"
        )
        obs = create_feature_table(annotations, level="class")
        pd.testing.assert_frame_equal(exp, obs)

    def test_create_feature_table_sub_class(self):
        exp = pd.DataFrame(
            {
                "ARSENIC": [0, 0, 0, 0, 1],
                "BETA-LACTAM": [1, 0, 1, 1, 0],
                "CEPHALOSPORIN": [0, 1, 0, 0, 0],
                "EFFLUX": [0, 0, 0, 0, 1],
            },
            index=["contig01", "contig02", "contig03", "contig08", "contig13"],
        )
        exp.index.name = "Contig id"
        exp.columns.name = "Subclass"
        annotations = AMRFinderPlusAnnotationsDirFmt(
            self.get_data_path("annotations_contigs_1"), mode="r"
        )
        obs = create_feature_table(annotations, level="subclass")
        pd.testing.assert_frame_equal(exp, obs)

    def test_key_error(self):
        annotations = AMRFinderPlusAnnotationsDirFmt(
            self.get_data_path("annotations_protein"), mode="r"
        )
        with self.assertRaisesRegex(KeyError, "solely from protein data"):
            create_feature_table(annotations)

    def test_empty_data_error(self):
        annotations = AMRFinderPlusAnnotationsDirFmt()
        with open(annotations.path / "sample1_amr_all_mutations.tsv", "w"):
            pass
        with self.assertRaisesRegex(ValueError, "File is empty"):
            create_feature_table(annotations)
