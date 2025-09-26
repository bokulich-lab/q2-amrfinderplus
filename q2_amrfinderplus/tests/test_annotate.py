from unittest.mock import MagicMock, call, patch

import qiime2
from q2_types.genome_data import GenesDirectoryFormat, ProteinsDirectoryFormat
from qiime2.plugin.testing import TestPluginBase

from q2_amrfinderplus.annotate import _annotate, annotate
from q2_amrfinderplus.types import (
    AMRFinderPlusAnnotationsDirFmt,
    AMRFinderPlusDatabaseDirFmt,
)


class TestAnnotate(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    def setUp(self):
        super().setUp()
        self.mags = qiime2.Artifact.import_data(
            "SampleData[MAGs]", self.get_data_path("sample_data_mags")
        )
        self.contigs = qiime2.Artifact.import_data(
            "SampleData[Contigs]", self.get_data_path("contigs")
        )
        self.proteins = qiime2.Artifact.import_data(
            "GenomeData[Proteins]", self.get_data_path("proteins_per_sample")
        )
        self.loci = qiime2.Artifact.import_data(
            "GenomeData[Loci]", self.get_data_path("loci_per_sample")
        )
        self.feature_data_mags = qiime2.Artifact.import_data(
            "FeatureData[MAG]", self.get_data_path("feature_data_mag")
        )

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
    def test__annotate(
        self,
        mock_create_empty_files,
        mock_run_amrfinderplus,
        mock_get_file_paths,
        mock_create_sample_dirs,
        mock_create_sample_dict,
    ):
        # Create mock for the AMRFinderPlusDatabaseDirFmt input
        amrfinderplus_db = AMRFinderPlusDatabaseDirFmt()

        # Call the function with mostly default inputs
        result = _annotate(amrfinderplus_db)

        # Ensure the output is the correct types
        self.assertIsInstance(result[0], AMRFinderPlusAnnotationsDirFmt)
        self.assertIsInstance(result[1], AMRFinderPlusAnnotationsDirFmt)
        self.assertIsInstance(result[2], GenesDirectoryFormat)
        self.assertIsInstance(result[3], ProteinsDirectoryFormat)

    @patch("q2_amrfinderplus.annotate._validate_inputs")
    def test_annotate_pipeline_mags(self, mock_validate_inputs):
        mock_action = MagicMock(
            side_effect=[
                lambda *args, **kwargs: (
                    "annotations",
                    "all_mutations",
                    "genes",
                    "proteins",
                ),
                lambda x: ("collated_annotations",),
                lambda x: ("collated_genes",),
                lambda x: ("collated_proteins",),
                lambda x, y: ({1: "partitioned_seqs"},),
                lambda x, y: ({1: "partitioned_proteins"},),
                lambda x, y: ({1: "partitioned_loci"},),
            ]
        )

        mock_ctx = MagicMock(get_action=mock_action)
        annotate(
            ctx=mock_ctx,
            sequences=self.mags,
            proteins=self.proteins,
            loci=self.contigs,
            amrfinderplus_db=AMRFinderPlusDatabaseDirFmt(),
        )
        assert mock_ctx.get_action.call_args_list == [
            call("amrfinderplus", "_annotate"),
            call("amrfinderplus", "collate_annotations"),
            call("types", "collate_genes"),
            call("types", "collate_proteins"),
            call("types", "partition_sample_data_mags"),
            call("types", "partition_proteins"),
            call("types", "partition_loci"),
        ]

    @patch("q2_amrfinderplus.annotate._validate_inputs")
    def test_annotate_pipeline_contigs_only(self, mock_validate_inputs):
        mock_action = MagicMock(
            side_effect=[
                lambda *args, **kwargs: (
                    "annotations",
                    "all_mutations",
                    "genes",
                    "proteins",
                ),
                lambda x: ("collated_annotations",),
                lambda x: ("collated_genes",),
                lambda x: ("collated_proteins",),
                lambda x, y: ({1: "partitioned_seqs"},),
            ]
        )

        mock_ctx = MagicMock(get_action=mock_action)
        annotate(
            ctx=mock_ctx,
            sequences=self.contigs,
            amrfinderplus_db=AMRFinderPlusDatabaseDirFmt(),
        )
        assert mock_ctx.get_action.call_args_list == [
            call("amrfinderplus", "_annotate"),
            call("amrfinderplus", "collate_annotations"),
            call("types", "collate_genes"),
            call("types", "collate_proteins"),
            call("assembly", "partition_contigs"),
        ]

    @patch("q2_amrfinderplus.annotate._validate_inputs")
    def test_annotate_pipeline_feature_data_mags(self, mock_validate_inputs):
        mock_action = MagicMock(
            side_effect=[
                lambda *args, **kwargs: (
                    "annotations",
                    "all_mutations",
                    "genes",
                    "proteins",
                ),
                lambda x: ("collated_annotations",),
                lambda x: ("collated_genes",),
                lambda x: ("collated_proteins",),
                lambda x, y: ({1: "partitioned_seqs"},),
            ]
        )

        mock_ctx = MagicMock(get_action=mock_action)
        annotate(
            ctx=mock_ctx,
            sequences=self.feature_data_mags,
            amrfinderplus_db=AMRFinderPlusDatabaseDirFmt(),
            organism="Escherichia",
        )
        assert mock_ctx.get_action.call_args_list == [
            call("amrfinderplus", "_annotate"),
            call("amrfinderplus", "collate_annotations"),
            call("types", "collate_genes"),
            call("types", "collate_proteins"),
            call("types", "partition_feature_data_mags"),
        ]

    @patch("q2_amrfinderplus.annotate._validate_inputs")
    def test_annotate_pipeline_protein_only(self, mock_validate_inputs):
        mock_action = MagicMock(
            side_effect=[
                lambda *args, **kwargs: (
                    "annotations",
                    "all_mutations",
                    "genes",
                    "proteins",
                ),
                lambda x: ("collated_annotations",),
                lambda x: ("collated_genes",),
                lambda x: ("collated_proteins",),
                lambda x, y: ({1: "partitioned_seqs"},),
            ]
        )

        mock_ctx = MagicMock(get_action=mock_action)
        annotate(
            ctx=mock_ctx,
            proteins=self.proteins,
            amrfinderplus_db=AMRFinderPlusDatabaseDirFmt(),
        )
        assert mock_ctx.get_action.call_args_list == [
            call("amrfinderplus", "_annotate"),
            call("amrfinderplus", "collate_annotations"),
            call("types", "collate_genes"),
            call("types", "collate_proteins"),
            call("types", "partition_proteins"),
        ]
