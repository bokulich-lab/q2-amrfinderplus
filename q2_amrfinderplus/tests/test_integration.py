from pathlib import Path

import pandas as pd
import qiime2
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.genome_data import (
    GenesDirectoryFormat,
    LociDirectoryFormat,
    ProteinsDirectoryFormat,
)
from q2_types.per_sample_sequences import ContigSequencesDirFmt, MultiMAGSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import amrfinderplus
from skbio.io import read

from q2_amrfinderplus.types import (
    AMRFinderPlusAnnotationsDirFmt,
    AMRFinderPlusDatabaseDirFmt,
)


class TestAnnotatePipelineIntegration(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    def setUp(self):
        super().setUp()
        self.amrfinderplus_db = qiime2.Artifact.import_data(
            "AMRFinderPlusDatabase",
            AMRFinderPlusDatabaseDirFmt(self.get_data_path("minimal_database"), "r"),
        )
        self.contigs = qiime2.Artifact.import_data(
            "SampleData[Contigs]",
            ContigSequencesDirFmt(self.get_data_path("minimal_contigs"), "r"),
        )
        self.proteins = qiime2.Artifact.import_data(
            "GenomeData[Proteins]",
            ProteinsDirectoryFormat(self.get_data_path("minimal_proteins"), "r"),
        )
        self.organism_proteins = qiime2.Artifact.import_data(
            "GenomeData[Proteins]",
            ProteinsDirectoryFormat(
                self.get_data_path("minimal_organism_proteins"), "r"
            ),
        )
        self.loci = qiime2.Artifact.import_data(
            "GenomeData[Loci]",
            LociDirectoryFormat(self.get_data_path("minimal_loci"), "r"),
        )
        self.mags = qiime2.Artifact.import_data(
            "SampleData[MAGs]",
            MultiMAGSequencesDirFmt(
                self.get_data_path("minimal_sample_data_mags"), "r"
            ),
        )
        self.mag = qiime2.Artifact.import_data(
            "FeatureData[MAG]",
            MAGSequencesDirFmt(self.get_data_path("minimal_feature_data_mag"), "r"),
        )
        self.expected_stx_gene = self._read_fasta_sequences(
            self.get_data_path("minimal_contigs")
        )[0]
        self.expected_stx_gene_without_stop = self.expected_stx_gene[:-3]
        self.expected_stx_protein = self._read_fasta_sequences(
            self.get_data_path("minimal_proteins")
        )[0]
        self.expected_acrr_protein = self._read_fasta_sequences(
            self.get_data_path("minimal_organism_proteins")
        )[0]

    def test_annotate_contigs_proteins_loci(self):
        results = self._run_annotate_pipeline(
            sequences=self.contigs, proteins=self.proteins, loci=self.loci
        )

        hits = self._assert_stx_hits(results, method="EXACTP")
        self.assertEqual(set(hits["Protein id"]), {"contig1_1"})
        self.assertEqual(set(hits["Contig id"]), {"contig1"})
        self._assert_gene_sequences(results, self.expected_stx_gene)
        self._assert_protein_sequences(results, self.expected_stx_protein)

    def test_annotate_sample_data_mags(self):
        results = self._run_annotate_pipeline(sequences=self.mags)

        hits = self._assert_stx_hits(results, method="EXACTX")
        self.assertEqual(set(hits["Contig id"]), {"contig1"})
        self._assert_gene_sequences(results, self.expected_stx_gene_without_stop)

    def test_annotate_feature_data_mag(self):
        results = self._run_annotate_pipeline(sequences=self.mag)

        hits = self._assert_stx_hits(results, method="EXACTX")
        self.assertEqual(set(hits["Contig id"]), {"contig1"})
        self._assert_gene_sequences(results, self.expected_stx_gene_without_stop)

    def test_annotate_proteins_only(self):
        results = self._run_annotate_pipeline(proteins=self.proteins)

        hits = self._assert_stx_hits(results, method="EXACTP")
        self.assertEqual(set(hits["Protein id"]), {"contig1_1"})
        self._assert_protein_sequences(results, self.expected_stx_protein)

    def test_annotate_proteins_organism(self):
        results = self._run_annotate_pipeline(
            proteins=self.organism_proteins, organism="Escherichia"
        )

        annotations = self._annotation_frame(results.amr_annotations)
        self.assertEqual(len(annotations), 2)
        self.assertEqual(set(annotations["Protein id"]), {"acrR_mut"})
        self.assertEqual(set(annotations["Element symbol"]), {"acrR_T5N"})
        self.assertEqual(set(annotations["Method"]), {"POINTP"})
        self.assertEqual(
            set(annotations["Closest reference accession"]), {"WP_000101737.1"}
        )
        self._assert_protein_sequences(results, self.expected_acrr_protein)

        all_mutations = self._annotation_frame(results.amr_all_mutations)
        self.assertEqual(len(all_mutations), 6)
        self.assertEqual(set(all_mutations["Protein id"]), {"acrR_mut"})
        self.assertIn("acrR_T5N", set(all_mutations["Element symbol"]))
        self.assertEqual(set(all_mutations["Method"]), {"POINTP"})

    def test_annotate_contigs_parallel(self):
        with self.test_config:
            results = amrfinderplus.pipelines.annotate.parallel(
                self.amrfinderplus_db,
                sequences=self.contigs,
                proteins=self.proteins,
                loci=self.loci,
                plus=True,
                threads=1,
                num_partitions=2,
            )._result()
        self._assert_valid_results(results)

        hits = self._assert_stx_hits(results, method="EXACTP")
        self.assertEqual(set(hits["Protein id"]), {"contig1_1"})
        self.assertEqual(set(hits["Contig id"]), {"contig1"})
        self._assert_gene_sequences(results, self.expected_stx_gene)
        self._assert_protein_sequences(results, self.expected_stx_protein)

    def _run_annotate_pipeline(self, **kwargs):
        results = amrfinderplus.pipelines.annotate(
            amrfinderplus_db=self.amrfinderplus_db,
            plus=True,
            threads=1,
            num_partitions=1,
            **kwargs,
        )
        self._assert_valid_results(results)
        return results

    def _assert_valid_results(self, results):
        expected_formats = {
            "amr_annotations": AMRFinderPlusAnnotationsDirFmt,
            "amr_all_mutations": AMRFinderPlusAnnotationsDirFmt,
            "amr_genes": GenesDirectoryFormat,
            "amr_proteins": ProteinsDirectoryFormat,
        }
        for name, expected_format in expected_formats.items():
            result = getattr(results, name)
            result.validate()
            self.assertIs(result.format, expected_format)

    def _assert_stx_hits(self, results, method):
        hits = self._annotation_frame(results.amr_annotations)

        self.assertEqual(len(hits), 2)
        self.assertEqual(set(hits["Element symbol"]), {"stxA2"})
        self.assertEqual(set(hits["Method"]), {method})
        self.assertEqual(set(hits["Closest reference accession"]), {"AAA16360.1"})
        return hits

    def _assert_gene_sequences(self, results, expected):
        genes = results.amr_genes.view(GenesDirectoryFormat)
        self.assertEqual(self._read_fasta_sequences(genes.path), [expected, expected])

    def _assert_protein_sequences(self, results, expected):
        proteins = results.amr_proteins.view(ProteinsDirectoryFormat)
        self.assertEqual(
            self._read_fasta_sequences(proteins.path), [expected, expected]
        )

    def _annotation_frame(self, artifact):
        annotations = artifact.view(AMRFinderPlusAnnotationsDirFmt)
        frames = [
            pd.read_csv(path, sep="\t")
            for path in annotations.path.rglob("*_amr_*.tsv")
        ]
        return pd.concat(frames, ignore_index=True)

    def _read_fasta_sequences(self, directory):
        sequences = []
        for path in sorted(Path(directory).rglob("*.fasta")):
            sequences.extend(
                str(sequence) for sequence in read(str(path), format="fasta")
            )
        return sequences
