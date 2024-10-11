# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tempfile

import pandas as pd
import qiime2
from pandas._testing import assert_frame_equal
from qiime2.core.exceptions import ValidationError
from qiime2.plugin.testing import TestPluginBase

from q2_amrfinderplus.types._format import (
    AMRFinderPlusAnnotationFormat,
    AMRFinderPlusAnnotationsDirFmt,
    AMRFinderPlusDatabaseDirFmt,
    _create_path,
)
from q2_amrfinderplus.types._transformer import (
    _metadata_transformer_helper,
    combine_dataframes,
)


class TestAMRFinderPlusTypesAndFormats(TestPluginBase):
    package = "q2_amrfinderplus.types.tests"

    def test_amrfinderplus_database_directory_format_validate_positive(self):
        format = AMRFinderPlusDatabaseDirFmt(self.get_data_path("database"), mode="r")
        format.validate()

    def test_amrfinderplus_annotation_format_validate_positive(self):
        filepath = self.get_data_path(
            "annotation/no_coordinates/"
            "aa447c99-ecd9-4c4a-a53b-4df6999815dd_amr_annotations.tsv"
        )

        format = AMRFinderPlusAnnotationFormat(filepath, mode="r")
        format.validate()

    def test_amrfinderplus_annotation_format_validate_positive_coordinates(self):
        filepath = self.get_data_path(
            "annotation/coordinates/"
            "e026af61-d911-4de3-a957-7e8bf837f30d_amr_annotations.tsv"
        )
        format = AMRFinderPlusAnnotationFormat(filepath, mode="r")
        format.validate()

    def test_amrfinderplus_annotation_format_validate_positive_empty(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_file_path = os.path.join(temp_dir, "amr_annotations.tsv")
            with open(temp_file_path, "w"):
                pass
            format = AMRFinderPlusAnnotationFormat(temp_file_path, mode="r")
            format.validate()

    def test_amrfinderplus_annotation_format_validation_error(self):
        with self.assertRaises(ValidationError) as context:
            path = self.get_data_path("annotation_wrong/amr_annotation.tsv")
            format = AMRFinderPlusAnnotationFormat(path, mode="r")
            format.validate()

            header_coordinates = [
                "Protein identifier",
                "Contig id",
                "Start",
                "Stop",
                "Strand",
                "Gene symbol",
                "Sequence name",
                "Scope",
                "Element type",
                "Element subtype",
                "Class",
                "Subclass",
                "Method",
                "Target length",
                "Reference sequence length",
                "% Coverage of reference sequence",
                "% Identity to reference sequence",
                "Alignment length",
                "Accession of closest sequence",
                "Name of closest sequence",
                "HMM id",
                "HMM description",
            ]
            expected_message = (
                "Header line does not match AMRFinderPlusAnnotation format. Must "
                "consist of the following values: "
                + ", ".join(header_coordinates)
                + ".\nWhile Contig id, Start, Stop and Strand are optional."
                + "\n\nFound instead: "
                + "Incorrect Header 1, Incorrect Header 2, Incorrect Header 3"
            )

            self.assertEqual(str(context.exception), expected_message)

    def test_amrfinderplus_annotations_dir_fmt_feature(self):
        dirpath = self.get_data_path("annotation/coordinates")
        annotations = AMRFinderPlusAnnotationsDirFmt(dirpath, mode="r")
        assert isinstance(annotations, AMRFinderPlusAnnotationsDirFmt)

    def test_amrfinderplus_annotations_dir_fmt_sample(self):
        dirpath = self.get_data_path("annotation")
        annotations = AMRFinderPlusAnnotationsDirFmt(dirpath, mode="r")
        assert isinstance(annotations, AMRFinderPlusAnnotationsDirFmt)

    def test_amrfinderplus_annotations_dir_fmt_path_maker_dir_name(self):
        fmt = AMRFinderPlusAnnotationsDirFmt()
        path = fmt.annotations_path_maker(
            name="annotations", id="id", dir_name="dir_name"
        )
        self.assertEqual(
            str(path), os.path.join(str(fmt), "dir_name/id_amr_annotations.tsv")
        )

    def test_amrfinderplus_annotations_dir_fmt_path_maker(self):
        fmt = AMRFinderPlusAnnotationsDirFmt()
        path = fmt.annotations_path_maker(name="annotations", id="id")
        self.assertEqual(str(path), os.path.join(str(fmt), "id_amr_annotations.tsv"))

    def test__create_path_annotations_absolute(self):
        dir_format = AMRFinderPlusAnnotationsDirFmt(
            self.get_data_path("annotation/coordinates"), "r"
        )
        file_path = (
            dir_format.path / "e026af61-d911-4de3-a957-7e8bf837f30d_amr_annotations.tsv"
        )
        path, _id = _create_path(path=file_path, relative=False, dir_format=dir_format)

        self.assertEqual(path, str(file_path))
        self.assertEqual(_id, "e026af61-d911-4de3-a957-7e8bf837f30d")

    def test__create_path_all_mutations_relative(self):
        dir_format = AMRFinderPlusAnnotationsDirFmt(
            self.get_data_path("all_mutations/coordinates"), "r"
        )
        file_path = (
            dir_format.path
            / "e026af61-d911-4de3-a957-7e8bf837f30d_amr_all_mutations.tsv"
        )
        path, _id = _create_path(path=file_path, relative=True, dir_format=dir_format)

        self.assertEqual(
            path, "e026af61-d911-4de3-a957-7e8bf837f30d_amr_all_mutations.tsv"
        )
        self.assertEqual(_id, "e026af61-d911-4de3-a957-7e8bf837f30d")

    def test_amrfinderplus_annotation_dirfmt_samples_annotation_dict(self):
        annotations = AMRFinderPlusAnnotationsDirFmt(
            self.get_data_path("annotation"), mode="r"
        )

        obs = annotations.annotation_dict()
        exp = {
            "coordinates": {
                "e026af61-d911-4de3-a957-7e8bf837f30d": os.path.join(
                    str(annotations),
                    "coordinates/"
                    "e026af61-d911-4de3-a957-7e8bf837f30d_amr_annotations.tsv",
                ),
            },
            "no_coordinates": {
                "aa447c99-ecd9-4c4a-a53b-4df6999815dd": os.path.join(
                    str(annotations),
                    "no_coordinates/"
                    "aa447c99-ecd9-4c4a-a53b-4df6999815dd_amr_annotations.tsv",
                ),
            },
        }
        self.assertDictEqual(obs, exp)

    def test_amrfinderplus_annotation_dirfmt_annotation_dict(self):
        annotations = AMRFinderPlusAnnotationsDirFmt(
            self.get_data_path("annotation/coordinates"), mode="r"
        )

        obs = annotations.annotation_dict()
        exp = {
            "e026af61-d911-4de3-a957-7e8bf837f30d": os.path.join(
                str(annotations),
                "e026af61-d911-4de3-a957-7e8bf837f30d_amr_annotations.tsv",
            )
        }

        self.assertDictEqual(obs, exp)


class MetadataTransformerUtilsTest(TestPluginBase):
    package = "q2_amrfinderplus.types.tests"

    def setUp(self):
        super().setUp()
        # Setup test data
        self.df1 = pd.DataFrame(
            {
                "Sample/MAG_ID": ["id_value_1", "id_value_2"],
                "col1": ["val1", "val3"],
                "col2": ["val2", "val4"],
            }
        )

        self.df2 = pd.DataFrame(
            {
                "Sample/MAG_ID": ["id_value_3", "id_value_4"],
                "col1": ["val5", "val7"],
                "col2": ["val6", "val8"],
            }
        )

        self.expected_combined_df = pd.DataFrame(
            {
                "Sample/MAG_ID": [
                    "id_value_1",
                    "id_value_2",
                    "id_value_3",
                    "id_value_4",
                ],
                "col1": ["val1", "val3", "val5", "val7"],
                "col2": ["val2", "val4", "val6", "val8"],
            }
        )

        self.expected_combined_df.index = self.expected_combined_df.index.astype(str)
        self.expected_combined_df.index.name = "id"

    def test_combine_dataframes(self):
        # Test combine_dataframes function
        df_list = [self.df1, self.df2]
        combined_df = combine_dataframes(df_list)
        pd.testing.assert_frame_equal(combined_df, self.expected_combined_df)


class TestAMRFinderPlusTransformers(TestPluginBase):
    package = "q2_amrfinderplus.types.tests"

    def test_annotations_feature_data_mags_transformer_helper(self):
        self._test_helper("annotations_feature_data_mags", "feature_data.tsv")

    def test_annotations_sample_data_contigs_transformer_helper(self):
        self._test_helper("annotations_sample_data_contigs", "sample_data_contigs.tsv")

    def test_annotations_sample_data_mags_transformer_helper(self):
        self._test_helper("annotations_sample_data_mags", "sample_data_mags.tsv")

    def test_mutations_feature_data_mags_transformer_helper(self):
        self._test_helper("mutations_feature_data_mags", "feature_data.tsv")

    def test_mutations_sample_data_contigs_transformer_helper(self):
        self._test_helper("mutations_sample_data_contigs", "sample_data_contigs.tsv")

    def test_mutations_sample_data_mags_transformer_helper(self):
        self._test_helper("mutations_sample_data_mags", "sample_data_mags.tsv")

    def _test_helper(self, data, table_name):
        df_expected = pd.read_csv(
            self.get_data_path(f"metadata_tables/{table_name}"),
            sep="\t",
        )
        df_expected.index = df_expected.index.astype(str)
        df_expected.index.name = "id"
        dir_fmt = AMRFinderPlusAnnotationsDirFmt(self.get_data_path(data), "r")
        df_obs = _metadata_transformer_helper(dir_fmt)
        assert_frame_equal(df_expected, df_obs)

    def test_annotations_sample_data_mags_to_Metadata(self):
        transformer = self.get_transformer(
            AMRFinderPlusAnnotationsDirFmt, qiime2.Metadata
        )
        fmt = AMRFinderPlusAnnotationsDirFmt(
            self.get_data_path("annotations_sample_data_mags"), "r"
        )

        metadata_obt = transformer(fmt)

        self.assertIsInstance(metadata_obt, qiime2.Metadata)
