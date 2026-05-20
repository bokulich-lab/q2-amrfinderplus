# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
from collections import defaultdict

import pandas as pd
from q2_types.feature_data import MixedCaseDNAFASTAFormat, ProteinFASTAFormat
from qiime2.core.exceptions import ValidationError
from qiime2.plugin import model

HEADER_SHARED_COLUMNS = {
    "Scope",
    "Class",
    "Subclass",
    "Method",
    "Target length",
    "Reference sequence length",
    "Alignment length",
    "HMM description",
    "Hierarchy node",
}
GENE_SYMBOL_COLUMNS = {"Gene symbol", "Element symbol"}


class TextFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass


class BinaryFormat(model.BinaryFileFormat):
    def _validate_(self, level):
        pass


class AMRFinderPlusDatabaseDirFmt(model.DirectoryFormat):
    amr_lib = model.File("AMR.LIB", format=TextFormat)
    amr_lib_comp = model.FileCollection(r"^AMR\.LIB\.h3.$", format=BinaryFormat)
    amrprot = model.File("AMRProt.fa", format=ProteinFASTAFormat)
    amrprot_blast = model.FileCollection(r"^AMRProt\.fa\.p..$", format=BinaryFormat)
    amrprot_mutation = model.File("AMRProt-mutation.tsv", format=TextFormat)
    amrprot_suppress = model.File("AMRProt-suppress.tsv", format=TextFormat)
    amrprot_susceptible = model.File("AMRProt-susceptible.tsv", format=TextFormat)
    amrprot_susceptible_fasta = model.File(
        "AMRProt-susceptible.fa", format=ProteinFASTAFormat
    )
    amr_cds = model.File("AMR_CDS.fa", format=MixedCaseDNAFASTAFormat)
    amr_cds_blast = model.FileCollection(r"^AMR_CDS\.fa\.n..$", format=BinaryFormat)
    fam = model.File("fam.tsv", format=TextFormat)
    taxgroup = model.File("taxgroup.tsv", format=TextFormat)
    version = model.File("version.txt", format=TextFormat)
    db_fmt_version = model.File("database_format_version.txt", format=TextFormat)
    amr_dna = model.FileCollection(
        r"^AMR_DNA-[a-zA-Z_]+\.fa$", format=MixedCaseDNAFASTAFormat
    )
    amr_dna_comp = model.FileCollection(
        r"^AMR_DNA-[a-zA-Z_]+\.fa\.n..$", format=BinaryFormat
    )
    amr_dna_tsv = model.FileCollection(r"^AMR_DNA-[a-zA-Z_]+\.tsv$", format=TextFormat)

    @amr_lib_comp.set_path_maker
    def amr_lib_comp_path_maker(self, extension):
        return "AMR.LIB.%s" % extension

    @amrprot_blast.set_path_maker
    def amrprot_blast_path_maker(self, extension):
        return "AMRProt.%s" % extension

    @amr_cds_blast.set_path_maker
    def amr_cds_blast_path_maker(self, extension):
        return "AMR_CDS.%s" % extension

    @amr_dna.set_path_maker
    def amr_dna_path_maker(self, species):
        return "AMR_DNA-%s" % species

    @amr_dna_comp.set_path_maker
    def amr_dna_comp_path_maker(self, species, extension):
        return "AMR_DNA-%s.%s" % (species, extension)

    @amr_dna_tsv.set_path_maker
    def amr_dna_tsv_path_maker(self, species):
        return "AMR_DNA-%s.tsv" % species


class AMRFinderPlusAnnotationFormat(model.TextFileFormat):
    def _validate(self):
        try:
            header_obs = pd.read_csv(str(self), sep="\t", nrows=0).columns.tolist()
            header_obs_set = set(header_obs)
            missing_columns = HEADER_SHARED_COLUMNS - header_obs_set
            has_gene_symbol_column = not GENE_SYMBOL_COLUMNS.isdisjoint(header_obs_set)

            if missing_columns or not has_gene_symbol_column:
                raise ValidationError(
                    "Header line does not match AMRFinderPlusAnnotationFormat. Must "
                    "contain the required output columns and at least one of: "
                    + ", ".join(sorted(GENE_SYMBOL_COLUMNS))
                    + "."
                    + "\n\nMissing columns: "
                    + ", ".join(sorted(missing_columns))
                    + "\n\nFound instead: "
                    + ", ".join(header_obs)
                )
        except pd.errors.EmptyDataError:
            pass

    def _validate_(self, level):
        self._validate()


class AMRFinderPlusAnnotationsDirFmt(model.DirectoryFormat):
    annotations = model.FileCollection(
        r".*amr_(annotations|all_mutations)\.tsv$", format=AMRFinderPlusAnnotationFormat
    )

    def annotation_dict(self, relative=False):
        """
        For per sample directories it returns a mapping of sample id to
        another dictionary where keys represent the file name and values
        correspond to the filepath for each file.
        For files, it returns a mapping of file name to filepath for each file.
        The suffixes "_amr_annotations" and "_amr_all_mutations" are removed from
        filenames.

        Parameters
        ---------
        relative : bool
            Whether to return filepaths relative to the directory's location.
            Returns absolute filepaths by default.

        Returns
        -------
        dict
            Mapping of filename -> filepath as described above.
            Or mapping of sample id -> dict {filename: filepath} as
            described above.
            Both levels of the dictionary are sorted alphabetically by key.
        """
        ids = defaultdict(dict)
        for entry in self.path.iterdir():
            if entry.is_dir():
                outer_id = entry.name
                for path in entry.iterdir():
                    file_path, inner_id = _create_path(
                        path=path, relative=relative, dir_format=self
                    )

                    ids[outer_id][inner_id] = str(file_path)
                ids[outer_id] = dict(sorted(ids[outer_id].items()))
            else:
                file_path, inner_id = _create_path(
                    path=entry, relative=relative, dir_format=self
                )

                ids[inner_id] = str(file_path)

        return dict(sorted(ids.items()))

    @annotations.set_path_maker
    def annotations_path_maker(self, name, id, dir_name=""):
        return os.path.join(dir_name, f"{id}_amr_{name}.tsv")


def _create_path(path, relative, dir_format):
    """
    This function processes the input file path to generate an absolute or relative
    path string and the sample or MAG ID derived from the file name. The ID is
    extracted by removing the suffix "_amr_annotations" or "_amr_all_mutations" from the
    file name. The created path and ID are used to build the annotation_dict that maps
    IDs to filepaths.

    Parameters:
    ---------
        path : Path
            A Path object representing the file path to process.
        relative : bool
            A flag indicating whether the returned path should be relative
            to the directory formats path or absolute.
        dir_format : AMRFinderplusAnnotationDirFmt
            An object of class "AMRFinderplusAnnotationDirFmt".

    Returns:
    -------
        path_dict : str
            The full relative or absolut path to the file.
        _id : str
            The sample or MAG ID derived from the file name.
    """
    file_name = path.stem

    # Remove suffix from filename to create id
    if file_name.endswith("_amr_annotations"):
        _id = file_name[:-16]
    else:
        _id = file_name[:-18]

    path_dict = (
        path.absolute().relative_to(dir_format.path.absolute())
        if relative
        else path.absolute()
    )
    return str(path_dict), _id
