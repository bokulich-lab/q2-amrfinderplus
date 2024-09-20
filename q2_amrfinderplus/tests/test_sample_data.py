import os
from unittest.mock import MagicMock, patch, mock_open

import pandas as pd
from q2_types.genome_data import ProteinsDirectoryFormat, LociDirectoryFormat
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase

from q2_amrfinderplus.sample_data import annotate_sample_data_amrfinderplus, \
    create_count_table, read_in_txt
from q2_amrfinderplus.types import AMRFinderPlusDatabaseDirFmt


def mock_run_amrfinderplus_n(
        amrfinderplus_db,
        dna_sequences,
        protein_sequences,
        gff,
        organism,
        plus,
        report_all_equal,
        ident_min,
        curated_ident,
        coverage_min,
        translation_table,
        annotation_format,
        report_common,
        threads,
        amr_annotations_path,
        amr_genes_path,
        amr_proteins_path,
        amr_all_mutations_path
):
    with open(amr_annotations_path, "w"):
        pass
    if organism:
        with open(amr_all_mutations_path, "w"):
            pass
    if dna_sequences:
        with open(amr_genes_path, "w"):
            pass
    if protein_sequences:
        with open(amr_proteins_path, "w"):
            pass


class TestAnnotateSampleDataAMRFinderPlus(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    def test_annotate_sample_data_amrfinderplus_mags(self):
        mags = MultiMAGSequencesDirFmt()
        os.mkdir(mags.path / "sample_1")
        with open(mags.path / "sample_1" / "mag_1.fasta", "w"):
            pass
        self._helper(mags=mags, organism="Escherichia")

    def test_annotate_sample_data_amrfinderplus_mags_proteins_loci(self):
        mags = MultiMAGSequencesDirFmt()
        proteins = ProteinsDirectoryFormat()
        loci = LociDirectoryFormat()
        os.mkdir(mags.path / "sample_1")
        with open(mags.path / "sample_1" / "mag_1.fasta", "w"):
            pass
        os.mkdir(proteins.path / "sample_1")
        with open(proteins.path / "sample_1" / "mag_1.fasta", "w"):
            pass
        os.mkdir(loci.path / "sample_1")
        with open(loci.path / "sample_1" / "mag_1.gff", "w"):
            pass

        self._helper(mags=mags, proteins=proteins, loci=loci)

    def test_annotate_sample_data_amrfinderplus_proteins(self):
        proteins = ProteinsDirectoryFormat()
        os.mkdir(proteins.path / "sample_1")
        with open(proteins.path / "sample_1" / "mag_1.fasta", "w"):
            pass
        self._helper(proteins=proteins)

    def _helper(self, mags=None, organism=None, proteins=None, loci=None):
        amrfinderplus_db = AMRFinderPlusDatabaseDirFmt()
        mock_create_count_table = MagicMock()
        mock_read_in_txt = MagicMock()
        with patch(
                "q2_amrfinderplus.sample_data.run_amrfinderplus_analyse",
                side_effect=mock_run_amrfinderplus_n,
        ), patch(
            "q2_amrfinderplus.sample_data.read_in_txt", mock_read_in_txt
        ), patch(
            "q2_amrfinderplus.sample_data.create_count_table",
            mock_create_count_table,
        ):
            (amr_annotations, amr_all_mutations, amr_genes,
             amr_proteins, feature_table) = annotate_sample_data_amrfinderplus(
                mags=mags,
                proteins=proteins,
                loci=loci,
                amrfinderplus_db=amrfinderplus_db,
                organism=organism,
            )
            self.assertTrue(
                os.path.exists(amr_annotations.path / "sample_1" /
                               "mag_1_amr_annotations.tsv")
            )

            if mags:
                self.assertTrue(
                    os.path.exists(amr_genes.path / "sample_1" /
                                   "mag_1_amr_genes.fasta")
                )
            else:
                self.assertTrue(os.path.exists(amr_genes.path / "empty.fasta"))
            if organism:
                self.assertTrue(
                    os.path.exists(amr_all_mutations.path / "sample_1" /
                                   "mag_1_amr_all_mutations.tsv")
                )
            else:
                self.assertTrue(
                    os.path.exists(
                        amr_all_mutations.path / "empty_amr_all_mutations.tsv")
                )
            if proteins:
                self.assertTrue(
                    os.path.exists(amr_proteins.path / "sample_1" /
                                   "mag_1_amr_proteins.fasta")
                )
            else:
                self.assertTrue(os.path.exists(amr_proteins.path / "empty.fasta"))


class TestReadInTxt(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    @patch("builtins.open", new_callable=mock_open,
           read_data="col1\tcol2\nA\t1\nB\t2\nA\t3")
    @patch("pandas.read_csv")
    def test_read_in_txt(self, mock_read_csv, mock_file):
        # Mock data that would be read by pd.read_csv
        mock_df = pd.DataFrame({
            'col1': ['A', 'B', 'A']
        })

        # Mock the behavior of pd.read_csv
        mock_read_csv.return_value = mock_df

        # Test parameters
        path = "dummy_path.txt"
        sample_mag_id = "Sample123"
        column_name = "col1"

        # Call the function under test
        result = read_in_txt(path, sample_mag_id, column_name)

        # Expected DataFrame after processing
        expected_df = pd.DataFrame({
            column_name: ['A', 'B'],
            sample_mag_id: ['2', '1']  # value_counts result for A is 2, for B is 1
        })

        # Convert to string to match the function's output format
        expected_df = expected_df.astype(str)

        # Check that the file was read using pandas
        mock_read_csv.assert_called_once_with(path, sep="\t")

        # Assert that the returned DataFrame matches the expected output
        pd.testing.assert_frame_equal(result, expected_df)


class TestCreateCountTable(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    @classmethod
    def setUpClass(cls):
        cls.gene_count_df = [pd.DataFrame(
            {
                "ARO Term": ["mdtF", "mgrA", "OprN", "mepA"],
                "sample1": ["1", "1", "1", "1"],
            }
        ),
            pd.DataFrame(
                {
                    "ARO Term": ["mdtE", "mgrA", "OprN", "mepA"],
                    "sample2": ["1", "1", "1", "1"],
                }
            )
            ]

        cls.frequency_table = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2"],
                "OprN": ["1", "1"],
                "mdtE": ["0", "1"],
                "mdtF": ["1", "0"],
                "mepA": ["1", "1"],
                "mgrA": ["1", "1"],

            }
        )
        cls.frequency_table.set_index("sample_id", inplace=True)

    def test_create_count_table(self):
        # Create observed count table with create_count_table function
        obs = create_count_table(self.gene_count_df)
        obs = obs.astype(str)

        # Define expected count table
        exp = self.frequency_table

        # Compare expected and observed count table
        pd.testing.assert_frame_equal(exp, obs)

    def test_create_count_table_value_error(self):
        # Assert if ValueError is called when empy list is passed
        self.assertRaises(ValueError, create_count_table, [])