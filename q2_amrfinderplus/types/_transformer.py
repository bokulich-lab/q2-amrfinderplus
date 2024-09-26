# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
from pathlib import Path

import pandas as pd
import qiime2

from q2_amrfinderplus.plugin_setup import plugin
from q2_amrfinderplus.types import AMRFinderPlusAnnotationsDirFmt


@plugin.register_transformer
def _1(data: AMRFinderPlusAnnotationsDirFmt) -> qiime2.Metadata:
    return qiime2.Metadata(_transfomer_helper(data))


def _transfomer_helper(data):
    df_list = []
    for file_dir_name in os.listdir(str(data)):
        # Check the directory structure
        if os.path.isdir(os.path.join(str(data), file_dir_name)):
            for file in glob.glob(os.path.join(str(data), file_dir_name, "*")):
                file_name = Path(file).stem

                # Annotations file from sample data mags
                if file_name.endswith("_amr_annotations"):
                    id_value = file_dir_name + "/" + file_name[:-16]

                # Mutations file from sample data mags
                else:
                    id_value = file_dir_name + "/" + file_name[:-18]

                create_append_df(
                    file_path=file,
                    df_list=df_list,
                    id_value=id_value,
                )
        else:
            # Annotations file from feature data mags or sample data contigs
            if file_dir_name.endswith("_amr_annotations.tsv"):
                id_value = file_dir_name[:-20]
            # Mutations file from feature data mags
            else:
                id_value = file_dir_name[:-22]

            create_append_df(
                file_path=os.path.join(str(data), file_dir_name),
                df_list=df_list,
                id_value=id_value,
            )

    return combine_dataframes(df_list)


def create_append_df(file_path, df_list, id_value):
    # Read in df
    df = pd.read_csv(file_path, sep="\t")

    # Insert column with sample or mag IDs
    df.insert(0, "Sample/MAG_ID", id_value)

    # Append df to df list
    df_list.append(df)


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
