import pandas as pd

from q2_amrfinderplus.types import AMRFinderPlusAnnotationsDirFmt


def create_feature_table(
    annotations: AMRFinderPlusAnnotationsDirFmt,
) -> pd.DataFrame:
    df = pd.DataFrame()
    sample_dict = annotations.annotation_dict()

    if type(next(iter(sample_dict.values()), None)) == str:
        sample_dict = {"": sample_dict}

    # Iterate over files in directory
    for sample_id, file_dict in sample_dict.items():
        for _id, file_fp in file_dict.items():
            try:
                file_df = pd.read_csv(
                    filepath_or_buffer=file_fp,
                    sep="\t",
                    usecols=["Contig id", "Gene symbol", "Start", "Stop", "Strand"],
                )
            except ValueError():
                raise (ValueError(f"Error reading file: {file_fp}"))

            df = pd.concat([df, file_df])

    df.drop_duplicates(keep="first", inplace=True)
    df_pivot = pd.crosstab(df["Contig id"], df["Gene symbol"])

    return df_pivot
