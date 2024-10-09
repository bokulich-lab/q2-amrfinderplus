# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_types.genome_data import GenomeData
from qiime2.core.type import SemanticType

AMRFinderPlusDatabase = SemanticType("AMRFinderPlusDatabase")
AMRFinderPlusAnnotations = SemanticType(
    "AMRFinderPlusAnnotations",
    variant_of=GenomeData.field["type"],
)
