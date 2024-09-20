import os
from functools import reduce

import pandas as pd
from q2_types.genome_data import (
    GenesDirectoryFormat,
    LociDirectoryFormat,
    ProteinsDirectoryFormat,
)
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt

from q2_amrfinderplus.types import (
    AMRFinderPlusAnnotationsDirFmt,
    AMRFinderPlusDatabaseDirFmt,
)
from q2_amrfinderplus.utils import _validate_inputs, run_amrfinderplus_analyse


def annotate_sample_data_amrfinderplus(
    amrfinderplus_db: AMRFinderPlusDatabaseDirFmt,
    mags: MultiMAGSequencesDirFmt = None,
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
    pd.DataFrame,
):
    # Validate input and parameter combinations
    _validate_inputs(
        mags, loci, proteins, ident_min, curated_ident, report_common, plus, organism
    )

    # Set up common parameters for run_amrfinderplus_analyse
    common_params = locals().copy()
    del common_params["mags"]
    del common_params["proteins"]
    del common_params["loci"]

    # Innit output formats
    amr_annotations = AMRFinderPlusAnnotationsDirFmt()
    amr_all_mutations = AMRFinderPlusAnnotationsDirFmt()
    amr_genes = GenesDirectoryFormat()
    amr_proteins = ProteinsDirectoryFormat()
    frequency_list = []

    # Create iterator for samples with sample_dict
    if mags:
        sample_iterator = mags.sample_dict().items()
    else:
        # Monkey patch the sample_dict instance method of MultiMAGSequencesDirFmt to
        # ProteinsDirectoryFormat because it should have the same per sample structure
        proteins.pathspec = r".+\.(fa|faa|fasta)$"
        proteins.sample_dict = MultiMAGSequencesDirFmt.sample_dict.__get__(
            proteins, ProteinsDirectoryFormat
        )
        sample_iterator = proteins.sample_dict().items()

    # Iterate over paths of MAGs
    for sample_id, files_dict in sample_iterator:
        # Create sample directories in output directories
        os.mkdir(f"{amr_annotations}/{sample_id}")
        if mags:
            os.mkdir(f"{amr_genes}/{sample_id}")
        if proteins:
            os.mkdir(f"{amr_proteins}/{sample_id}")
        if organism:
            os.mkdir(f"{amr_all_mutations}/{sample_id}")

        for mag_id, file_fp in files_dict.items():
            # Run amrfinderplus
            run_amrfinderplus_analyse(
                # dna_sequences path is the mag full path if mags are specified and
                # None if no mags are specifies
                dna_sequences=file_fp if mags else None,
                # protein_sequences path is constructed if mags and proteins are
                # specified, the mag full path is used when only proteins is specified.
                # If only mags are specified and not proteins, the path is None.
                protein_sequences=proteins.path / sample_id / f"{mag_id}.fasta"
                if mags and proteins
                else file_fp
                if not mags
                else None,
                gff=loci.path / sample_id / f"{mag_id}.gff" if loci else None,
                amr_annotations_path=amr_annotations.path
                / sample_id
                / f"{mag_id}_amr_annotations.tsv",
                amr_genes_path=amr_genes.path / sample_id / f"{mag_id}_amr_genes.fasta",
                amr_proteins_path=amr_proteins.path
                / sample_id
                / f"{mag_id}_amr_proteins.fasta",
                amr_all_mutations_path=amr_all_mutations.path
                / sample_id
                / f"{mag_id}_amr_all_mutations.tsv",
                **common_params,
            )

            # Create frequency dataframe and append it to list
            frequency_df = read_in_txt(
                path=amr_annotations.path / sample_id / f"{mag_id}_amr_annotations.tsv",
                sample_mag_id=f"{sample_id}/{mag_id}",
                column_name="Gene symbol",
            )
            frequency_list.append(frequency_df)

    feature_table = create_count_table(df_list=frequency_list)

    if not mags:
        with open(os.path.join(str(amr_genes), "empty.fasta"), "w"):
            pass

    if not proteins:
        with open(os.path.join(str(amr_proteins), "empty.fasta"), "w"):
            pass

    if not organism:
        with open(
            os.path.join(str(amr_all_mutations), "empty_amr_all_mutations.tsv"), "w"
        ):
            pass

    return (amr_annotations, amr_all_mutations, amr_genes, amr_proteins, feature_table)


def read_in_txt(path: str, sample_mag_id: str, column_name: str):
    # Read in txt file to pd.Dataframe
    df = pd.read_csv(path, sep="\t")

    # Process the df, create count table
    df = df[column_name].value_counts().reset_index()
    df.columns = [column_name, sample_mag_id]

    df = df.astype(str)
    return df


def create_count_table(df_list: list) -> pd.DataFrame:
    # Remove all empty lists from df_list
    df_list = [df for df in df_list if not df.empty]

    # Raise ValueError if df_list is empty. This happens when no ARGs were detected
    if not df_list:
        raise ValueError(
            "No AMR genes could be identified and no output can be created."
        )

    # Merge all dfs contained in df_list
    df = reduce(
        lambda left, right: pd.merge(left, right, on=left.columns[0], how="outer"),
        df_list,
    )

    # Process the df to meet all requirements for a FeatureTable
    df = df.transpose()
    df = df.fillna(0)
    df.columns = df.iloc[0]
    df = df.drop(df.index[0])
    df.columns.name = None
    df.index.name = "sample_id"
    return df
