import os
import subprocess

from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import ContigSequencesDirFmt, MultiMAGSequencesDirFmt
from qiime2.util import duplicate

from q2_amrfinderplus.types import AMRFinderPlusAnnotationsDirFmt

EXTERNAL_CMD_WARNING = (
    "Running external command line application(s). "
    "This may print messages to stdout and/or stderr.\n"
    "The command(s) being run are below. These commands "
    "cannot be manually re-run as they will depend on "
    "temporary files that no longer exist."
)


def run_command(cmd, cwd=None, verbose=True):
    if verbose:
        print(EXTERNAL_CMD_WARNING)
        print("\nCommand:", end=" ")
        print(" ".join(cmd), end="\n\n")
    subprocess.run(cmd, check=True, cwd=cwd, stderr=subprocess.PIPE)


def _validate_inputs(
    sequences, loci, proteins, ident_min, curated_ident, report_common, plus, organism
):
    # Ensure that at least sequences or proteins is provided
    if not sequences and not proteins:
        raise ValueError('"--i-sequences" or "--i-proteins" input has to be provided.')

    # Check if loci is provided with sequences but without proteins
    # (invalid combination)
    if sequences and loci and not proteins:
        raise ValueError(
            '"--i-loci" input can only be given in combination with "--i-proteins" '
            "input."
        )

    # Check if sequences and proteins are provided together but without loci
    # (invalid combination)
    if sequences and not loci and proteins:
        raise ValueError(
            '"--i-sequences" and "--i-proteins" inputs together can only '
            'be given in combination with "--i-loci" input.'
        )

    # Validate that ident_min and curated_ident are not used together
    if ident_min and curated_ident:
        raise ValueError(
            '"--p-ident-min" and "--p-curated-ident" cannot be used simultaneously.'
        )

    # Check that report_common is only used with plus and organism
    if report_common and (not plus or not organism):
        raise ValueError('"--p-report-common" requires "--p-plus" and "--p-organism".')


def _run_amrfinderplus_analyse(
    amrfinderplus_db,
    dna_path,
    protein_path,
    gff_path,
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
    amr_genes_path=None,
    amr_proteins_path=None,
    amr_all_mutations_path=None,
):
    cmd = [
        "amrfinder",
        "--database",
        str(amrfinderplus_db),
        "-o",
        amr_annotations_path,
        "--print_node",
    ]
    # Creates nucleotide fasta output if DNA sequences are given as input
    if dna_path:
        cmd.extend(
            [
                "-n",
                dna_path,
                "--nucleotide_output",
                amr_genes_path,
            ]
        )
    # Creates protein fasta output if protein sequences are given as input
    if protein_path:
        cmd.extend(
            [
                "-p",
                protein_path,
                "--protein_output",
                amr_proteins_path,
            ]
        )
    if gff_path:
        cmd.extend(["-g", gff_path])
    if threads:
        cmd.extend(["--threads", str(threads)])
    # Creates all mutations output if an organism is specified
    if organism:
        cmd.extend(
            [
                "--organism",
                organism,
                "--mutation_all",
                amr_all_mutations_path,
            ]
        )
    if plus:
        cmd.append("--plus")
    if report_all_equal:
        cmd.append("--report_all_equal")
    if ident_min:
        cmd.extend(["--ident_min", str(ident_min)])
    if curated_ident:
        cmd.extend(["--ident_min", "-1"])
    if coverage_min:
        cmd.extend(["--coverage_min", str(coverage_min)])
    if translation_table:
        cmd.extend(["--translation_table", str(translation_table)])
    if annotation_format:
        cmd.extend(["--annotation_format", str(annotation_format)])
    if report_common:
        cmd.append("--report_common")
    if organism in [
        "Acinetobacter",
        "Burkholderia_cepacia_complex",
        "Escherichia_coli_Shigella",
        "Klebsiella",
        "Serratia",
    ]:
        cmd.append("--gpipe_org")

    try:
        run_command(cmd=cmd)
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode("utf-8"))
        if "gff_check.cpp" in e.stderr.decode("utf-8"):
            raise Exception(
                "GFF file error: Either there is data missing in one GFF file or an "
                "incorrect GFF input format was specified with the parameter "
                "'--p-annotation-format'. Please check https://github.com/ncbi/amr/wiki"
                "/Running-AMRFinderPlus#input-file-formats for documentation and "
                "choose the correct GFF input format. Please inspect stdout and stderr "
                "to learn more."
            )
        else:
            raise Exception(
                "An error was encountered while running AMRFinderPlus, "
                f"(return code {e.returncode}), please inspect stdout and stderr to "
                "learn more."
            )


def _create_empty_files(
    sequences, proteins, organism, amr_genes, amr_proteins, amr_all_mutations
):
    # Creates empty files in output artifacts amr_genes, amr_proteins and
    # amr_all_mutations because artifacts can not be empty
    if not sequences:
        with open(os.path.join(str(amr_genes), "empty.fasta"), "w"):
            pass
        print(
            colorify(
                '"amr_genes" output is empty because no "--i-sequences" input '
                "was given."
            )
        )

    if not proteins:
        with open(os.path.join(str(amr_proteins), "empty.fasta"), "w"):
            pass
        print(
            colorify(
                '"amr_proteins" output is empty because no "--i-proteins" input '
                "was given."
            )
        )

    if not organism:
        with open(
            os.path.join(str(amr_all_mutations), "empty_amr_all_mutations.tsv"), "w"
        ):
            pass
        print(
            colorify(
                '"amr_all_mutations" output is empty because no "--p-organism" '
                "parameter was given."
            )
        )


def _create_sample_dirs(
    sequences,
    proteins,
    organism,
    amr_annotations,
    amr_genes,
    amr_proteins,
    amr_all_mutations,
    sample_id,
):
    os.makedirs(os.path.join(str(amr_annotations), sample_id), exist_ok=True)
    if sequences:
        os.makedirs(os.path.join(str(amr_genes), sample_id), exist_ok=True)
    if proteins:
        os.makedirs(os.path.join(str(amr_proteins), sample_id), exist_ok=True)
    if organism:
        os.makedirs(os.path.join(str(amr_all_mutations), sample_id), exist_ok=True)


def _create_sample_dict(proteins, sequences):
    if sequences:
        # For SampleData[MAGs]
        if isinstance(sequences, MultiMAGSequencesDirFmt):
            sample_dict = sequences.sample_dict()

        # For SampleData[Contigs]
        elif isinstance(sequences, ContigSequencesDirFmt):
            file_dict = sequences.sample_dict()
            # Create fake sample for sample_dict
            sample_dict = {"": file_dict}

        # For FeatureData[MAG]
        elif isinstance(sequences, MAGSequencesDirFmt):
            file_dict = sequences.feature_dict()
            # Create fake sample for sample_dict
            sample_dict = {"": file_dict}

    else:
        # For GenomeData[Proteins] with per sample directory structure
        if any(item.is_dir() for item in proteins.path.iterdir()):
            sample_dict = proteins.file_dict()

        # For GenomeData[Proteins] with no per sample directory structure
        else:
            file_dict = proteins.file_dict()
            # Create sample_dict with fake sample
            sample_dict = {"": file_dict}

    return sample_dict


def _get_file_paths(sequences, proteins, loci, _id, file_fp, sample_id=""):
    # If mags is provided
    if sequences:
        dna_path = file_fp

        # If proteins are provided, construct the expected protein file path.
        if proteins:
            protein_path = os.path.join(str(proteins), sample_id, f"{_id}.fasta")

            # Raise an error if the expected protein file does not exist.
            if not os.path.exists(protein_path):
                raise ValueError(
                    f"Proteins file for ID '{_id}' is missing in proteins input."
                )
        else:
            protein_path = None

    # If only proteins are provided (without mags), determine dna and protein file path.
    else:
        dna_path = None
        protein_path = file_fp

    # If loci are provided, construct the expected GFF file path.
    if loci:
        gff_path = os.path.join(str(loci), sample_id, f"{_id}.gff")

        # Raise an error if the expected GFF file does not exist.
        if not os.path.exists(gff_path):
            raise ValueError(f"GFF file for ID '{_id}' is missing in loci input.")
    else:
        gff_path = None

    return dna_path, protein_path, gff_path


def colorify(string: str):
    return "%s%s%s" % ("\033[1;33m", string, "\033[0m")


def collate_annotations(
    annotations: AMRFinderPlusAnnotationsDirFmt,
) -> AMRFinderPlusAnnotationsDirFmt:
    collated_annotations = AMRFinderPlusAnnotationsDirFmt()

    for annotation in annotations:
        for item in annotation.path.iterdir():
            target = collated_annotations.path / item.name
            if item.is_file():
                duplicate(item, target)
            elif item.is_dir():
                target.mkdir(exist_ok=True)
                for file in item.iterdir():
                    duplicate(file, target / file.name)
    return collated_annotations
