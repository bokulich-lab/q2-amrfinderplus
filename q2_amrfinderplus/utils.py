import subprocess


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


def _validate_inputs(mags, loci, proteins, ident_min, curated_ident, report_common,
                     plus, organism):
    # Check if loci is provided with mags but without proteins (invalid combination)
    if mags and loci and not proteins:
        raise ValueError(
            '"--i-loci" input can only be given in combination with "--i-proteins" '
            'input.'
        )

    # Check if mags and proteins are provided together but without loci
    # (invalid combination)
    if mags and not loci and proteins:
        raise ValueError(
            '"--i-mags" and "--i-proteins" inputs together can only '
            'be given in combination with "--i-loci" input.'
        )

    # Ensure that at least mags or proteins is provided
    if not mags and not proteins:
        raise ValueError('"--i-mags" or "--i-proteins" input has to be provided.')

    # Validate that ident_min and curated_ident are not used together
    if ident_min and curated_ident:
        raise ValueError('"--p-ident-min" and "--p-curated-ident" cannot be used '
                         'simultaneously.')

    # Check that report_common is only used with plus and organism
    if report_common and (not plus or not organism):
        raise ValueError('"--p-report-common" requires "--p-plus" and "--p-organism".')


def run_amrfinderplus_analyse(
        amrfinderplus_db,
        dna_sequences,
        protein_sequences,
        gff,
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
        amr_all_mutations_path=None
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
    if dna_sequences:
        cmd.extend(
            [
                "-n",
                dna_sequences,
                "--nucleotide_output",
                amr_genes_path,
            ]
        )
    # Creates protein fasta output if protein sequences are given as input
    if protein_sequences:
        cmd.extend(
            [
                "-p",
                protein_sequences,
                "--protein_output",
                amr_proteins_path,
            ]
        )
    if gff:
        cmd.extend(["-g", gff])
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
    if organism in ["Acinetobacter",
                    "Burkholderia_cepacia_complex",
                    "Escherichia_coli_Shigella",
                    "Klebsiella",
                    "Serratia"]:
        cmd.append("--gpipe_org")

    try:
        run_command(cmd=cmd)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running AMRFinderPlus, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
