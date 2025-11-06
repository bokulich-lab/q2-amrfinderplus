import os
from typing import Union

from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.genome_data import (
    GenesDirectoryFormat,
    LociDirectoryFormat,
    ProteinsDirectoryFormat,
)
from q2_types.per_sample_sequences import ContigSequencesDirFmt, MultiMAGSequencesDirFmt

from q2_amrfinderplus.types import (
    AMRFinderPlusAnnotationsDirFmt,
    AMRFinderPlusDatabaseDirFmt,
)
from q2_amrfinderplus.utils import (
    _create_empty_files,
    _create_sample_dict,
    _create_sample_dirs,
    _get_file_paths,
    _run_amrfinderplus_analyse,
    _validate_inputs,
    remove_duplicate_ids_fasta,
)


def annotate(
    amrfinderplus_db: AMRFinderPlusDatabaseDirFmt,
    sequences: Union[
        MultiMAGSequencesDirFmt, ContigSequencesDirFmt, MAGSequencesDirFmt
    ] = None,
    proteins: ProteinsDirectoryFormat = None,
    loci: LociDirectoryFormat = None,
    organism: str = None,
    plus: bool = False,
    report_all_equal: bool = False,
    ident_min: float = None,
    curated_ident: bool = False,
    coverage_min: float = 0.5,
    translation_table: str = "11",
    annotation_format: str = "prodigal",
    report_common: bool = False,
    threads: int = None,
) -> (
    AMRFinderPlusAnnotationsDirFmt,
    AMRFinderPlusAnnotationsDirFmt,
    GenesDirectoryFormat,
    ProteinsDirectoryFormat,
):
    # Validate input and parameter combinations
    _validate_inputs(
        sequences,
        loci,
        proteins,
        ident_min,
        curated_ident,
        report_common,
        plus,
        organism,
    )

    # Set up common parameters for _run_amrfinderplus_analyse
    common_params = {
        k: v for k, v in locals().items() if k not in ("sequences", "proteins", "loci")
    }

    # Innit output formats
    amr_annotations = AMRFinderPlusAnnotationsDirFmt()
    amr_all_mutations = AMRFinderPlusAnnotationsDirFmt()
    amr_genes = GenesDirectoryFormat()
    amr_proteins = ProteinsDirectoryFormat()

    # Create sample_dict to iterate over input files
    sample_dict = _create_sample_dict(proteins, sequences)

    # Iterate over sample_dict
    for sample_id, files_dict in sample_dict.items():
        # Create sample directories in output directories
        _create_sample_dirs(
            sequences,
            proteins,
            organism,
            amr_annotations,
            amr_genes,
            amr_proteins,
            amr_all_mutations,
            sample_id,
        )

        for _id, file_fp in files_dict.items():
            # Construct and validate file input paths for amrfinderplus
            dna_path, protein_path, gff_path = _get_file_paths(
                sequences,
                proteins,
                loci,
                _id,
                file_fp,
                sample_id,
            )

            # Define paths for output files
            amr_annotations_path = os.path.join(
                str(amr_annotations), sample_id, f"{_id}_amr_annotations.tsv"
            )
            amr_genes_path = os.path.join(
                str(amr_genes), sample_id, f"{_id}_amr_genes.fasta"
            )
            amr_proteins_path = os.path.join(
                str(amr_proteins), sample_id, f"{_id}_amr_proteins.fasta"
            )
            amr_all_mutations_path = os.path.join(
                str(amr_all_mutations), sample_id, f"{_id}_amr_all_mutations.tsv"
            )

            # Run amrfinderplus
            _run_amrfinderplus_analyse(
                dna_path=dna_path,
                protein_path=protein_path,
                gff_path=gff_path,
                amr_annotations_path=amr_annotations_path,
                amr_genes_path=amr_genes_path,
                amr_proteins_path=amr_proteins_path,
                amr_all_mutations_path=amr_all_mutations_path,
                **common_params,
            )

            # Remove duplicate sequence IDs from gene and protein fasta outputs if
            # report-all-equal is set to True
            if report_all_equal and proteins:
                remove_duplicate_ids_fasta(amr_proteins_path)
            if report_all_equal and sequences:
                remove_duplicate_ids_fasta(amr_genes_path)

    # Create empty files for empty output artifacts if needed
    _create_empty_files(
        sequences, proteins, organism, amr_genes, amr_proteins, amr_all_mutations
    )

    return amr_annotations, amr_all_mutations, amr_genes, amr_proteins
