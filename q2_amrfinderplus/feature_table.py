import pandas as pd

from q2_amrfinderplus.types import AMRFinderPlusAnnotationsDirFmt


def create_feature_table(
    annotations: AMRFinderPlusAnnotationsDirFmt,
    level: str = "gene",
) -> pd.DataFrame:
    df = pd.DataFrame()
    sample_dict = annotations.annotation_dict()

    level_mapping = {
        "class": "Class",
        "subclass": "Subclass",
    }

    # Check if sample_dict is nested and create fake sample if needed
    if type(next(iter(sample_dict.values()), None)) == str:
        sample_dict = {"": sample_dict}

    # Loop over all files, read in dataframes and concatenate them
    for sample_id, file_dict in sample_dict.items():
        for _id, file_fp in file_dict.items():
            try:
                if level == "gene":
                    level_column = _get_gene_symbol_column(file_fp)
                else:
                    level_column = level_mapping[level]

                file_df = pd.read_csv(
                    filepath_or_buffer=file_fp,
                    sep="\t",
                    usecols=[
                        "Contig id",
                        level_column,
                        "Start",
                        "Stop",
                        "Strand",
                    ],
                )
            except pd.errors.EmptyDataError as e:
                raise ValueError(
                    "File is empty. All mutations output is empty if no organism was "
                    f"specified.\n\nOriginal error: {e}"
                )
            except ValueError as e:
                raise ValueError(
                    "If the annotations were created solely from protein data, there "
                    "is no positional information and no gene abundance per contig "
                    f"can be calculated.\n\nOriginal error: {e}"
                )

            df = pd.concat([df, file_df])

    # Drop duplicated rows and pivot table
    df.drop_duplicates(keep="first", inplace=True)
    df_pivot = pd.crosstab(df["Contig id"], df[level_column])

    return df_pivot


def _get_gene_symbol_column(file_fp):
    header = pd.read_csv(file_fp, sep="\t", nrows=0).columns
    if "Gene symbol" in header:
        return "Gene symbol"
    return "Element symbol"
