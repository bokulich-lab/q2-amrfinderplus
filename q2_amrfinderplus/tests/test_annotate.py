from unittest.mock import patch

from q2_types.genome_data import GenesDirectoryFormat, ProteinsDirectoryFormat
from qiime2.plugin.testing import TestPluginBase

from q2_amrfinderplus.annotate import annotate
from q2_amrfinderplus.types import (
    AMRFinderPlusAnnotationsDirFmt,
    AMRFinderPlusDatabaseDirFmt,
)


class TestAnnotate(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    @patch("q2_amrfinderplus.annotate._validate_inputs")
    @patch(
        "q2_amrfinderplus.annotate._create_sample_dict",
        return_value={"sample1": {"id1": "file_path"}},
    )
    @patch("q2_amrfinderplus.annotate._create_sample_dirs")
    @patch(
        "q2_amrfinderplus.annotate._get_file_paths",
        return_value=("dna_path", "protein_path", "gff_path"),
    )
    @patch("q2_amrfinderplus.annotate._run_amrfinderplus_analyse")
    @patch("q2_amrfinderplus.annotate._create_empty_files")
    def test_annotate(
        self,
        mock_create_empty_files,
        mock_run_amrfinderplus,
        mock_get_file_paths,
        mock_create_sample_dirs,
        mock_create_sample_dict,
        mock_validate_inputs,
    ):
        # Create mock for the AMRFinderPlusDatabaseDirFmt input
        amrfinderplus_db = AMRFinderPlusDatabaseDirFmt()

        # Call the function with mostly default inputs
        result = annotate(amrfinderplus_db)

        # Ensure the output is the correct types
        self.assertIsInstance(result[0], AMRFinderPlusAnnotationsDirFmt)
        self.assertIsInstance(result[1], AMRFinderPlusAnnotationsDirFmt)
        self.assertIsInstance(result[2], GenesDirectoryFormat)
        self.assertIsInstance(result[3], ProteinsDirectoryFormat)
