# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
import qiime2

from q2_amrfinderplus.plugin_setup import plugin
from q2_amrfinderplus.types import AMRFinderPlusAnnotationsDirFmt


@plugin.register_transformer
def _1(data: AMRFinderPlusAnnotationsDirFmt) -> qiime2.Metadata:
    return qiime2.Metadata(_metadata_transformer_helper(data))


def _metadata_transformer_helper(data):
    df_list = []

    if any(item.is_dir() for item in data.path.iterdir()):
        annotation_dict = data.annotation_dict()

    else:
        file_dict = data.annotation_dict()
        # Create annotation_dict with fake sample
        annotation_dict = {"": file_dict}

    for outer_id, files_dict in annotation_dict.items():
        for inner_id, file_fp in files_dict.items():
            id_value = f"{outer_id}/{inner_id}" if outer_id else inner_id

            # Create df and append to df_list
            df = pd.read_csv(file_fp, sep="\t")
            df.insert(0, "Sample/MAG_ID", id_value)
            df_list.append(df)

    return combine_dataframes(df_list)


def combine_dataframes(df_list):
    # Concat all dfs
    df_combined = pd.concat(df_list, axis=0)

    # Sort all values by sample/mag ID column
    df_combined.sort_values(by=df_combined.columns[0], inplace=True)

    # Reset and rename index and set it to string to conform to metadata format
    df_combined.reset_index(inplace=True, drop=True)
    df_combined.index.name = "id"
    df_combined.index = df_combined.index.astype(str)

    return df_combined
