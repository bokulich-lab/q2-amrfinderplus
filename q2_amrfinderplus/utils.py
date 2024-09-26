import os
import subprocess

from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.genome_data import ProteinsDirectoryFormat
from q2_types.per_sample_sequences import ContigSequencesDirFmt, MultiMAGSequencesDirFmt

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
    subprocess.run(cmd, check=True, cwd=cwd)


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
        str(amr_annotations_path),
        "--print_node",
    ]
    # Creates nucleotide fasta output if DNA sequences are given as input
    if dna_path:
        cmd.extend(
            [
                "-n",
                str(dna_path),
                "--nucleotide_output",
                str(amr_genes_path),
            ]
        )
    # Creates protein fasta output if protein sequences are given as input
    if protein_path:
        cmd.extend(
            [
                "-p",
                str(protein_path),
                "--protein_output",
                str(amr_proteins_path),
            ]
        )
    if gff_path:
        cmd.extend(["-g", str(gff_path)])
    if threads:
        cmd.extend(["--threads", str(threads)])
    # Creates all mutations output if an organism is specified
    if organism:
        cmd.extend(
            [
                "--organism",
                organism,
                "--mutation_all",
                str(amr_all_mutations_path),
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
        raise Exception(
            "An error was encountered while running AMRFinderPlus, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def _create_empty_files(
    sequences, proteins, organism, amr_genes, amr_proteins, amr_all_mutations
):
    # Creates empty files in output artifacts amr_genes, amr_proteins and
    # amr_all_mutations because artifacts can not be empty
    if not sequences:
        with open(amr_genes.path / "empty.fasta", "w"):
            pass
        print(
            colorify(
                '"amr_genes" output is empty because no "--i-sequences" input '
                "was given."
            )
        )

    if not proteins:
        with open(amr_proteins.path / "empty.fasta", "w"):
            pass
        print(
            colorify(
                '"amr_proteins" output is empty because no "--i-proteins" input '
                "was given."
            )
        )

    if not organism:
        with open(amr_all_mutations.path / "empty_amr_all_mutations.tsv", "w"):
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
    os.makedirs(amr_annotations.path / sample_id, exist_ok=True)
    if sequences:
        os.makedirs(amr_genes.path / sample_id, exist_ok=True)
    if proteins:
        os.makedirs(amr_proteins.path / sample_id, exist_ok=True)
    if organism:
        os.makedirs(amr_all_mutations.path / sample_id, exist_ok=True)


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
        proteins.pathspec = r".+\.(fa|faa|fasta)$"

        # Monkey patch the sample_dict instance method of MultiMAGSequencesDirFmt to
        # ProteinsDirectoryFormat if it has a sample data dir structure
        if any(item.is_dir() for item in proteins.path.iterdir()):
            proteins.sample_dict = MultiMAGSequencesDirFmt.sample_dict.__get__(
                proteins, ProteinsDirectoryFormat
            )
            sample_dict = proteins.sample_dict()
        # Monkey patch the feature_dict instance method of MAGSequencesDirFmt to
        # ProteinsDirectoryFormat if it has a feature data dir structure
        else:
            proteins.feature_dict = MAGSequencesDirFmt.feature_dict.__get__(
                proteins, ProteinsDirectoryFormat
            )
            file_dict = proteins.feature_dict()
            # create sample_dict with fake sample
            sample_dict = {"": file_dict}

    return sample_dict


def _get_file_paths(sequences, proteins, loci, id, file_fp, sample_id=""):
    # If mags is provided
    if sequences:
        dna_path = file_fp

        # If proteins are provided, construct the expected protein file path.
        if proteins:
            protein_path = proteins.path / sample_id / f"{id}.fasta"

            # Raise an error if the expected protein file does not exist.
            if not os.path.exists(protein_path):
                raise ValueError(
                    f"Proteins file for ID '{id}' is missing in proteins input."
                )
        else:
            protein_path = None

    # If only proteins are provided (without mags), determine dna and protein file path.
    else:
        dna_path = None
        protein_path = file_fp

    # If loci are provided, construct the expected GFF file path.
    if loci:
        gff_path = loci.path / sample_id / f"{id}.gff"

        # Raise an error if the expected GFF file does not exist.
        if not os.path.exists(gff_path):
            raise ValueError(f"GFF file for ID '{id}' is missing in loci input.")
    else:
        gff_path = None

    return dna_path, protein_path, gff_path


def colorify(string: str):
    return "%s%s%s" % ("\033[1;33m", string, "\033[0m")
