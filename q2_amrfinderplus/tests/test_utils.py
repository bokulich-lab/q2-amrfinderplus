import os
import subprocess
from pathlib import Path
from unittest.mock import MagicMock, call, patch

from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.genome_data import GenesDirectoryFormat, ProteinsDirectoryFormat
from q2_types.per_sample_sequences import ContigSequencesDirFmt, MultiMAGSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase

from q2_amrfinderplus.types import AMRFinderPlusAnnotationsDirFmt
from q2_amrfinderplus.utils import (
    EXTERNAL_CMD_WARNING,
    _create_empty_files,
    _create_sample_dict,
    _create_sample_dirs,
    _get_file_paths,
    _run_amrfinderplus_analyse,
    _validate_inputs,
    colorify,
    run_command,
)


class TestRunCommand(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    @patch("subprocess.run")
    @patch("builtins.print")
    def test_run_command_verbose(self, mock_print, mock_subprocess_run):
        # Mock command and working directory
        cmd = ["echo", "Hello"]
        cwd = "/test/directory"

        # Run the function with verbose=True
        run_command(cmd, cwd=cwd, verbose=True)

        # Check if subprocess.run was called with the correct arguments
        mock_subprocess_run.assert_called_once_with(cmd, check=True, cwd=cwd)

        # Check if the correct print statements were called
        mock_print.assert_has_calls(
            [
                call(EXTERNAL_CMD_WARNING),
                call("\nCommand:", end=" "),
                call("echo Hello", end="\n\n"),
            ]
        )

    @patch("subprocess.run")
    @patch("builtins.print")
    def test_run_command_non_verbose(self, mock_print, mock_subprocess_run):
        # Mock command and working directory
        cmd = ["echo", "Hello"]
        cwd = "/test/directory"

        # Run the function with verbose=False
        run_command(cmd, cwd=cwd, verbose=False)

        # Check if subprocess.run was called with the correct arguments
        mock_subprocess_run.assert_called_once_with(cmd, check=True, cwd=cwd)

        # Ensure no print statements were made
        mock_print.assert_not_called()


class TestRunAMRFinderPlusAnalyse(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    @patch("q2_amrfinderplus.utils.run_command")
    def test_run_amrfinderplus_analyse(self, mock_run_command):
        _run_amrfinderplus_analyse(
            amrfinderplus_db="amrfinderplus_db",
            dna_path="dna_sequences",
            protein_path="protein_sequences",
            gff_path="gff",
            organism="Acinetobacter",
            plus=True,
            report_all_equal=True,
            ident_min=1,
            curated_ident=False,
            coverage_min=1,
            translation_table="11",
            annotation_format="prodigal",
            report_common=True,
            threads=4,
            amr_annotations_path="amr_annotations_path",
            amr_genes_path="amr_genes_path",
            amr_proteins_path="amr_proteins_path",
            amr_all_mutations_path="amr_all_mutations_path",
        )
        mock_run_command.assert_called_once_with(
            cmd=[
                "amrfinder",
                "--database",
                "amrfinderplus_db",
                "-o",
                "amr_annotations_path",
                "--print_node",
                "-n",
                "dna_sequences",
                "--nucleotide_output",
                "amr_genes_path",
                "-p",
                "protein_sequences",
                "--protein_output",
                "amr_proteins_path",
                "-g",
                "gff",
                "--threads",
                "4",
                "--organism",
                "Acinetobacter",
                "--mutation_all",
                "amr_all_mutations_path",
                "--plus",
                "--report_all_equal",
                "--ident_min",
                "1",
                "--coverage_min",
                "1",
                "--translation_table",
                "11",
                "--annotation_format",
                "prodigal",
                "--report_common",
                "--gpipe_org",
            ],
        )

    @patch("q2_amrfinderplus.utils.run_command")
    def test_run_amrfinderplus_analyse_minimal(self, mock_run_command):
        _run_amrfinderplus_analyse(
            amrfinderplus_db="amrfinderplus_db",
            dna_path=None,
            protein_path=None,
            gff_path=None,
            organism=None,
            plus=False,
            report_all_equal=False,
            ident_min=None,
            curated_ident=True,
            coverage_min=None,
            translation_table=None,
            annotation_format=None,
            report_common=False,
            threads=None,
            amr_annotations_path="amr_annotations_path",
        )
        mock_run_command.assert_called_once_with(
            cmd=[
                "amrfinder",
                "--database",
                "amrfinderplus_db",
                "-o",
                "amr_annotations_path",
                "--print_node",
                "--ident_min",
                "-1",
            ],
        )

    @patch("q2_amrfinderplus.utils.run_command")
    def test_run_amrfinderplus_analyse_exception_message(self, mock_run_command):
        # Simulate subprocess.CalledProcessError
        mock_run_command.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd="amrfinder"
        )

        # Call the function and assert the exception message
        with self.assertRaises(Exception) as context:
            _run_amrfinderplus_analyse(
                amrfinderplus_db="mock_db",
                dna_path=None,
                protein_path=None,
                gff_path=None,
                organism=None,
                plus=False,
                report_all_equal=False,
                ident_min=None,
                curated_ident=False,
                coverage_min=0.5,
                translation_table="11",
                annotation_format="prodigal",
                report_common=False,
                threads=None,
                amr_annotations_path="mock_annotations_path",
            )

        # Assert the correct exception message is raised
        self.assertIn(
            "An error was encountered while running AMRFinderPlus",
            str(context.exception),
        )
        self.assertIn("(return code 1)", str(context.exception))


class TestValidateInputs(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    # Test when --i-loci is given without --i-proteins
    def test_loci_without_proteins(self):
        with self.assertRaisesRegex(
            ValueError, "can only be given in combination " 'with "--i-proteins"'
        ):
            _validate_inputs(
                sequences=True,
                loci=True,
                proteins=False,
                ident_min=None,
                curated_ident=None,
                report_common=None,
                plus=None,
                organism=None,
            )

    # Test when --i-mags and --i-proteins are given without --i-loci
    def test_mags_and_proteins_without_loci(self):
        with self.assertRaisesRegex(
            ValueError, "can only be given in combination " 'with "--i-loci"'
        ):
            _validate_inputs(
                sequences=True,
                loci=False,
                proteins=True,
                ident_min=None,
                curated_ident=None,
                report_common=None,
                plus=None,
                organism=None,
            )

    # Test when neither --i-mags nor --i-proteins is provided
    def test_missing_mags_and_proteins(self):
        with self.assertRaisesRegex(
            ValueError, '"--i-sequences" or "--i-proteins" input has to be provided'
        ):
            _validate_inputs(
                sequences=False,
                loci=False,
                proteins=False,
                ident_min=None,
                curated_ident=None,
                report_common=None,
                plus=None,
                organism=None,
            )

    # Test when both --p-ident-min and --p-curated-ident are given
    def test_ident_min_and_curated_ident(self):
        with self.assertRaisesRegex(
            ValueError,
            '"--p-ident-min" and '
            '"--p-curated-ident" cannot be used '
            "simultaneously",
        ):
            _validate_inputs(
                sequences=True,
                loci=None,
                proteins=None,
                ident_min=True,
                curated_ident=True,
                report_common=None,
                plus=None,
                organism=None,
            )

    # Test when --p-report-common is given but --p-plus or --p-organism is missing
    def test_report_common_without_plus_or_organism(self):
        with self.assertRaisesRegex(
            ValueError, '"--p-report-common" requires "--p-plus" and "--p-organism"'
        ):
            _validate_inputs(
                sequences=True,
                loci=None,
                proteins=None,
                ident_min=None,
                curated_ident=None,
                report_common=True,
                plus=False,
                organism=None,
            )


class TestGetFilePaths(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    @patch("os.path.exists")
    def test_mags_with_proteins_and_loci(self, mock_exists):
        # Mock the os.path.exists to simulate files existing
        mock_exists.side_effect = [True, True]  # First for protein, second for GFF

        # Call the function with mags, proteins, and loci
        dna_path, protein_path, gff_path = _get_file_paths(
            sequences=MagicMock(),
            proteins=MagicMock(path=Path("proteins")),
            loci=MagicMock(path=Path("loci")),
            _id="id",
            sample_id="sample1",
            file_fp="dna_file.fasta",
        )

        # Assertions
        self.assertEqual(dna_path, "dna_file.fasta")
        self.assertEqual(str(protein_path), "proteins/sample1/id.fasta")
        self.assertEqual(str(gff_path), "loci/sample1/id.gff")

    def test_mags_without_proteins_and_loci(self):
        # Call the function with mags, proteins, and loci
        dna_path, protein_path, gff_path = _get_file_paths(
            sequences=MagicMock(),
            proteins=None,
            loci=None,
            _id="sample123",
            file_fp="dna_file.fasta",
        )

        # Assertions
        self.assertEqual(dna_path, "dna_file.fasta")
        self.assertEqual(protein_path, None)
        self.assertEqual(gff_path, None)

    @patch("os.path.exists")
    def test_mags_with_missing_protein(self, mock_exists):
        # Mock os.path.exists to simulate the missing protein file
        mock_exists.side_effect = [False]  # Protein file does not exist

        # Call the function with mags and proteins, but no loci
        with self.assertRaisesRegex(
            ValueError, "Proteins file for ID 'sample123' is missing"
        ):
            _get_file_paths(
                sequences=MagicMock(),
                proteins=MagicMock(),
                loci=None,
                _id="sample123",
                sample_id="sample1",
                file_fp="dna_file.fasta",
            )

    @patch("os.path.exists")
    def test_loci_with_missing_gff(self, mock_exists):
        # Mock os.path.exists to simulate the protein file exists but GFF file is
        # missing
        mock_exists.side_effect = [False]  # Protein exists, GFF is missing

        # Call the function with proteins and loci, but no mags
        with self.assertRaises(ValueError) as context:
            _get_file_paths(
                sequences=None,
                proteins=None,
                loci=MagicMock(path=Path("/mock/loci/path")),
                _id="sample123",
                sample_id="sample1",
                file_fp="protein_file.fasta",
            )

        # Check that the exception message contains the correct text
        self.assertIn("GFF file for ID 'sample123' is missing", str(context.exception))


class TestCreateSampleDict(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    def test_create_sample_dict_sequences_multimags(self):
        dirpath = self.get_data_path("sample_data_mags")
        sequences = MultiMAGSequencesDirFmt(dirpath, "r")

        result = _create_sample_dict(proteins=None, sequences=sequences)

        self.assertEqual(
            result, {"sample1": {"mag": os.path.join(dirpath, "sample1", "mag.fasta")}}
        )

    def test_create_sample_dict_sequences_contigs(self):
        dirpath = self.get_data_path("contigs")
        sequences = ContigSequencesDirFmt(dirpath, "r")

        result = _create_sample_dict(proteins=None, sequences=sequences)

        self.assertEqual(
            result, {"": {"sample1": os.path.join(dirpath, "sample1_contigs.fasta")}}
        )

    def test_create_sample_dict_sequences_mag(self):
        dirpath = self.get_data_path("feature_data_mag")
        sequences = MAGSequencesDirFmt(dirpath, "r")

        result = _create_sample_dict(proteins=None, sequences=sequences)

        self.assertEqual(
            result,
            {
                "": {
                    "30ef72ed-84fd-4348-a418-9d68a9b88729": os.path.join(
                        dirpath, "30ef72ed-84fd-4348-a418-9d68a9b88729.fasta"
                    )
                }
            },
        )

    def test_create_sample_dict_proteins(self):
        dirpath = self.get_data_path("proteins")

        proteins = ProteinsDirectoryFormat(dirpath, "r")

        result = _create_sample_dict(proteins=proteins, sequences=None)

        self.assertEqual(
            result, {"": {"sample1": os.path.join(dirpath, "sample1.fasta")}}
        )

    def test_create_sample_dict_proteins_per_sample(self):
        dirpath = self.get_data_path("proteins_per_sample")

        proteins = ProteinsDirectoryFormat(dirpath, "r")

        result = _create_sample_dict(proteins=proteins, sequences=None)

        self.assertEqual(
            result,
            {"sample1": {"genome1": os.path.join(dirpath, "sample1", "genome1.fasta")}},
        )


class TestCreateEmptyFiles(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    def test_create_empty_files_all_false(self):
        amr_genes = GenesDirectoryFormat()
        amr_proteins = ProteinsDirectoryFormat()
        amr_all_mutations = AMRFinderPlusAnnotationsDirFmt()

        _create_empty_files(
            sequences=False,
            proteins=False,
            organism=False,
            amr_genes=amr_genes,
            amr_proteins=amr_proteins,
            amr_all_mutations=amr_all_mutations,
        )

        self.assertTrue(os.path.exists(os.path.join(str(amr_genes), "empty.fasta")))
        self.assertTrue(os.path.exists(os.path.join(str(amr_proteins), "empty.fasta")))
        self.assertTrue(
            os.path.exists(
                os.path.join(str(amr_all_mutations), "empty_amr_all_mutations.tsv")
            )
        )

    def test_create_empty_files_all_true(self):
        _create_empty_files(
            sequences=True,
            proteins=True,
            organism=True,
            amr_genes=GenesDirectoryFormat(),
            amr_proteins=ProteinsDirectoryFormat(),
            amr_all_mutations=AMRFinderPlusAnnotationsDirFmt(),
        )


class TestCreateSampleDirs(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    def test_create_sample_dirs_all_exist(self):
        amr_annotations = AMRFinderPlusAnnotationsDirFmt()
        amr_genes = GenesDirectoryFormat()
        amr_proteins = ProteinsDirectoryFormat()
        amr_all_mutations = AMRFinderPlusAnnotationsDirFmt()

        _create_sample_dirs(
            sequences=True,
            proteins=True,
            organism=True,
            amr_annotations=amr_annotations,
            amr_genes=amr_genes,
            amr_proteins=amr_proteins,
            amr_all_mutations=amr_all_mutations,
            sample_id="sample1",
        )

        self.assertTrue(os.path.exists(os.path.join(str(amr_annotations), "sample1")))
        self.assertTrue(os.path.exists(os.path.join(str(amr_genes), "sample1")))
        self.assertTrue(os.path.exists(os.path.join(str(amr_proteins), "sample1")))
        self.assertTrue(os.path.exists(os.path.join(str(amr_all_mutations), "sample1")))

    def test_create_sample_dirs_nothing(self):
        amr_annotations = AMRFinderPlusAnnotationsDirFmt()

        _create_sample_dirs(
            sequences=False,
            proteins=False,
            organism=False,
            amr_annotations=amr_annotations,
            amr_genes=None,
            amr_proteins=None,
            amr_all_mutations=None,
            sample_id="sample1",
        )

        self.assertTrue(os.path.exists(os.path.join(str(amr_annotations), "sample1")))


class TestColorify(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    def test_colorify(self):
        # Test if colorify wraps the string with the correct ANSI codes for yellow
        result = colorify("Hello")
        expected = "\033[1;33mHello\033[0m"
        self.assertEqual(result, expected)
