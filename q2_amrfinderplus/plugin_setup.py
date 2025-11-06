# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib

from q2_types.feature_data import FeatureData
from q2_types.feature_data_mag import MAG
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.genome_data import Genes, GenomeData, Loci, Proteins
from q2_types.per_sample_sequences import Contigs, MAGs
from q2_types.sample_data import SampleData
from qiime2.core.type import Bool, Choices, Float, Int, Range, Str
from qiime2.plugin import Citations, Plugin

from q2_amrfinderplus import __version__
from q2_amrfinderplus.annotate import annotate
from q2_amrfinderplus.database import fetch_amrfinderplus_db
from q2_amrfinderplus.feature_table import create_feature_table
from q2_amrfinderplus.types._format import (
    AMRFinderPlusAnnotationFormat,
    AMRFinderPlusAnnotationsDirFmt,
    AMRFinderPlusDatabaseDirFmt,
    BinaryFormat,
    TextFormat,
)
from q2_amrfinderplus.types._type import AMRFinderPlusAnnotations, AMRFinderPlusDatabase

citations = Citations.load("citations.bib", package="q2_amrfinderplus")

plugin = Plugin(
    name="amrfinderplus",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-amrfinderplus",
    package="q2_amrfinderplus",
    description=(
        "A plugin to find acquired antimicrobial resistance genes and point mutations "
        "in protein and/or assembled nucleotide sequences with NCBI-AMRFinderPlus."
    ),
    short_description="AMR annotation.",
    citations=[],
)

plugin.methods.register_function(
    function=fetch_amrfinderplus_db,
    inputs={},
    parameters={},
    outputs=[("amrfinderplus_db", AMRFinderPlusDatabase)],
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={
        "amrfinderplus_db": "AMRFinderPlus database.",
    },
    name="Download AMRFinderPlus database.",
    description="Download the latest version of the AMRFinderPlus database.",
    citations=[citations["feldgarden2021amrfinderplus"]],
)

organisms = [
    "Acinetobacter_baumannii",
    "Bordetella_pertussis",
    "Burkholderia_cepacia",
    "Burkholderia_mallei",
    "Burkholderia_pseudomallei",
    "Campylobacter",
    "Citrobacter_freundii",
    "Clostridioides_difficile",
    "Corynebacterium_diphtheriae",
    "Enterobacter_asburiae",
    "Enterobacter_cloacae",
    "Enterococcus_faecalis",
    "Enterococcus_faecium",
    "Escherichia",
    "Haemophilus_influenzae",
    "Klebsiella_oxytoca",
    "Klebsiella_pneumoniae",
    "Neisseria_gonorrhoeae",
    "Neisseria_meningitidis",
    "Pseudomonas_aeruginosa",
    "Salmonella",
    "Serratia_marcescens",
    "Staphylococcus_aureus",
    "Staphylococcus_pseudintermedius",
    "Streptococcus_agalactiae",
    "Streptococcus_pneumoniae",
    "Streptococcus_pyogenes",
    "Vibrio_cholerae",
    "Vibrio_parahaemolyticus",
    "Vibrio_vulnificus",
]


translation_tables = [
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "9",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "21",
    "22",
    "23",
    "24",
    "25",
    "26",
    "27",
    "28",
    "29",
    "30",
    "31",
    "33",
]


amrfinderplus_parameters = {
    "organism": Str % Choices(organisms),
    "plus": Bool,
    "report_all_equal": Bool,
    "ident_min": Float % Range(0, 1, inclusive_start=True, inclusive_end=True),
    "curated_ident": Bool,
    "coverage_min": Float % Range(0, 1, inclusive_start=True, inclusive_end=True),
    "translation_table": Str % Choices(translation_tables),
    "annotation_format": Str
    % Choices(
        "bakta",
        "genbank",
        "microscope",
        "patric",
        "pgap",
        "prodigal",
        "prokka",
        "pseudomonasdb",
        "rast",
        "standard",
    ),
    "report_common": Bool,
    "threads": Int % Range(0, None, inclusive_start=False),
}

amrfinderplus_parameter_descriptions = {
    "organism": (
        "Taxon used for screening known resistance causing point mutations and "
        "blacklisting of common, non-informative genes. Pathogen Detection taxgroup "
        "names can also be used."
    ),
    "plus": (
        "Provide results from 'Plus' genes such as virulence factors, stress-response "
        "genes, etc."
    ),
    "report_all_equal": (
        "Report all equally scoring BLAST and HMM matches. This will report multiple "
        "lines for a single element if there are multiple reference proteins that have "
        "the same score. On those lines the fields Accession of closest sequence and "
        "Name of closest sequence will be different showing each of the database "
        "proteins that are equally close to the query sequence."
    ),
    "ident_min": (
        "Minimum identity for a blast-based hit (Methods BLAST or PARTIAL). Setting "
        "this value to something other than -1 will override curated similarity "
        "cutoffs. We only recommend using this option if you have a specific reason."
    ),
    "curated_ident": (
        "Use the curated threshold for a blast-based hit, if it exists and 0.9 "
        "otherwise."
    ),
    "coverage_min": (
        "Minimum proportion of reference gene covered for a BLAST-based hit (Methods "
        "BLAST or PARTIAL)."
    ),
    "translation_table": "Translation table used for BLASTX.",
    "annotation_format": (
        "Specify the format of the GFF file in the loci input. 'standart' refers to "
        "NCBI resources such as GenBank and RefSeq."
    ),
    "report_common": "Report proteins common to a taxonomy group.",
    "threads": (
        "The number of threads to use for processing. AMRFinderPlus defaults to 4 on "
        "hosts with >= 4 cores. Setting this number higher than the number of cores on "
        "the running host may cause blastp to fail. Using more than 4 threads may "
        "speed up searches."
    ),
}

amrfinderplus_output_descriptions = {
    "amr_annotations": "Annotated AMR genes and mutations.",
    "amr_all_mutations": (
        "Report of genotypes at all locations screened for point mutations. These "
        "files allow you to distinguish between called point mutations that were the "
        "sensitive variant and the point mutations that could not be called because "
        "the sequence was not found. This file will contain all detected variants from "
        "the reference sequence, so it could be used as an initial screen for novel "
        "variants. Note 'Gene symbols' for mutations not in the database (identifiable "
        "by [UNKNOWN] in the Sequence name field) have offsets that are relative to "
        "the start of the sequence indicated in the field 'Accession of closest "
        "sequence' while 'Gene symbols' from known point-mutation sites have gene "
        "symbols that match the Pathogen Detection Reference Gene Catalog standardized "
        "nomenclature for point mutations."
    ),
    "amr_genes": (
        "Sequences that were identified by AMRFinderPlus as AMR genes. This will "
        "include the entire region that aligns to the references for point mutations."
    ),
    "amr_proteins": (
        "Protein Sequences that were identified by AMRFinderPlus as AMR genes. This "
        "will include the entire region that aligns to the references for point "
        "mutations"
    ),
}


amrfinderplus_input_descriptions = {
    "sequences": "MAGs or contigs to be annotated with AMRFinderPlus.",
    "proteins": "Protein sequences to be annotated with AMRFinderPlus.",
    "loci": (
        "GFF files to give sequence coordinates for proteins input. Required for "
        "combined searches of protein and DNA sequences."
    ),
    "amrfinderplus_db": "AMRFinderPlus Database.",
}


plugin.methods.register_function(
    function=annotate,
    inputs={
        "sequences": SampleData[MAGs | Contigs] | FeatureData[MAG],
        "proteins": GenomeData[Proteins],
        "loci": GenomeData[Loci],
        "amrfinderplus_db": AMRFinderPlusDatabase,
    },
    parameters=amrfinderplus_parameters,
    outputs=[
        ("amr_annotations", GenomeData[AMRFinderPlusAnnotations]),
        ("amr_all_mutations", GenomeData[AMRFinderPlusAnnotations]),
        ("amr_genes", GenomeData[Genes]),
        ("amr_proteins", GenomeData[Proteins]),
    ],
    input_descriptions=amrfinderplus_input_descriptions,
    parameter_descriptions=amrfinderplus_parameter_descriptions,
    output_descriptions=amrfinderplus_output_descriptions,
    name="Annotate MAGs or contigs with AMRFinderPlus.",
    description=(
        "Annotate MAGs or contigs with antimicrobial resistance genes with "
        "AMRFinderPlus. Check https://github.com/ncbi/amr/wiki for documentation."
    ),
    citations=[citations["feldgarden2021amrfinderplus"]],
)

plugin.methods.register_function(
    function=create_feature_table,
    inputs={"annotations": GenomeData[AMRFinderPlusAnnotations]},
    outputs=[("table", FeatureTable[Frequency])],
    parameters={},
    input_descriptions={"annotations": "AMR annotations."},
    output_descriptions={"table": "Frequency of AMR genes per contig."},
    parameter_descriptions={},
    name="Gene per contig frequency table",
    description=(
        "Create a gene per contig frequency table from AMRFinderPlus annotations."
    ),
)

plugin.register_semantic_type_to_format(
    AMRFinderPlusDatabase,
    artifact_format=AMRFinderPlusDatabaseDirFmt,
)
plugin.register_semantic_type_to_format(
    GenomeData[AMRFinderPlusAnnotations],
    artifact_format=AMRFinderPlusAnnotationsDirFmt,
)

plugin.register_formats(
    AMRFinderPlusDatabaseDirFmt,
    TextFormat,
    BinaryFormat,
    AMRFinderPlusAnnotationFormat,
    AMRFinderPlusAnnotationsDirFmt,
)

importlib.import_module("q2_amrfinderplus.types._transformer")
