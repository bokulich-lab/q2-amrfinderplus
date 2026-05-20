import pandas as pd

from q2_amrfinderplus.types import AMRFinderPlusAnnotationsDirFmt


def create_feature_table(
    annotations: AMRFinderPlusAnnotationsDirFmt,
    level: str = "gene",
) -> pd.DataFrame:
    frames = []
    sample_dict = annotations.annotation_dict()

    # Check if sample_dict is nested and create fake sample if needed
    if type(next(iter(sample_dict.values()), None)) == str:
        sample_dict = {"": sample_dict}

    # Loop over all files, read in dataframes and concatenate them
    for sample_id, file_dict in sample_dict.items():
        for _id, file_fp in file_dict.items():
            try:
                file_df = pd.read_csv(file_fp, sep="\t")

                if level == "gene":
                    if "Gene symbol" in file_df.columns:
                        level_column = "Gene symbol"
                    else:
                        level_column = "Element symbol"
                else:
                    level_column = level.capitalize()

                file_df = file_df[
                    ["Contig id", level_column, "Start", "Stop", "Strand"]
                ]

            except pd.errors.EmptyDataError as e:
                raise ValueError(
                    "File is empty. All mutations output is empty if no organism was "
                    f"specified.\n\nOriginal error: {e}"
                )
            except KeyError as e:
                raise KeyError(
                    "If the annotations were created solely from protein data, there "
                    "is no positional information and no gene abundance per contig "
                    f"can be calculated.\n\nOriginal error: {e}"
                )

            frames.append(file_df)

    # Drop duplicated rows and pivot table
    df = pd.concat(frames)
    df.drop_duplicates(keep="first", inplace=True)
    df_pivot = pd.crosstab(df["Contig id"], df[level_column])

    return df_pivot
