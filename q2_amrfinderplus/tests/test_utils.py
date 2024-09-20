from unittest.mock import patch

from qiime2.plugin.testing import TestPluginBase

from q2_amrfinderplus.utils import run_amrfinderplus_analyse, _validate_inputs


class TestRunAmrfinderplusAnalyse(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    @patch("q2_amrfinderplus.utils.run_command")
    def test_run_amrfinderplus_analyse(self, mock_run_command):
        run_amrfinderplus_analyse(
            amrfinderplus_db="amrfinderplus_db",
            dna_sequences="dna_sequences",
            protein_sequences="protein_sequences",
            gff="gff",
            organism="Escherichia",
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
            amr_all_mutations_path="amr_all_mutations_path"
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
                "Escherichia",
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
            ],
        )

    @patch("q2_amrfinderplus.utils.run_command")
    def test_run_amrfinderplus_analyse_minimal(self, mock_run_command):
        run_amrfinderplus_analyse(
            amrfinderplus_db="amrfinderplus_db",
            dna_sequences=None,
            protein_sequences=None,
            gff=None,
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


class TestValidateInputs(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    # Test when --i-loci is given without --i-proteins
    def test_loci_without_proteins(self):
        with self.assertRaisesRegex(ValueError, 'can only be given in combination '
                                                'with "--i-proteins"'):
            _validate_inputs(mags=True, loci=True, proteins=False, ident_min=None,
                             curated_ident=None, report_common=None, plus=None,
                             organism=None)

    # Test when --i-mags and --i-proteins are given without --i-loci
    def test_mags_and_proteins_without_loci(self):
        with self.assertRaisesRegex(ValueError, 'can only be given in combination '
                                                'with "--i-loci"'):
            _validate_inputs(mags=True, loci=False, proteins=True, ident_min=None,
                             curated_ident=None, report_common=None, plus=None,
                             organism=None)

    # Test when neither --i-mags nor --i-proteins is provided
    def test_missing_mags_and_proteins(self):
        with self.assertRaisesRegex(ValueError, '"--i-mags" or "--i-proteins" input '
                                                'has to be provided'):
            _validate_inputs(mags=False, loci=False, proteins=False, ident_min=None,
                             curated_ident=None, report_common=None, plus=None,
                             organism=None)

    # Test when both --p-ident-min and --p-curated-ident are given
    def test_ident_min_and_curated_ident(self):
        with self.assertRaisesRegex(ValueError, '"--p-ident-min" and '
                                                '"--p-curated-ident" cannot be used '
                                                'simultaneously'):
            _validate_inputs(mags=True, loci=None, proteins=None, ident_min=True,
                             curated_ident=True, report_common=None, plus=None,
                             organism=None)

    # Test when --p-report-common is given but --p-plus or --p-organism is missing
    def test_report_common_without_plus_or_organism(self):
        with self.assertRaisesRegex(ValueError, '"--p-report-common" requires '
                                                '"--p-plus" and "--p-organism"'):
            _validate_inputs(mags=True, loci=None, proteins=None, ident_min=None,
                             curated_ident=None, report_common=True, plus=False,
                             organism=None)
