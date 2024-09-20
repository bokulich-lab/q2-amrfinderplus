# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_types.sample_data import SampleData
from qiime2.plugin import Citations, Plugin

from q2_amrfinderplus import __version__
from q2_amrfinderplus.database import fetch_amrfinderplus_db
from q2_amrfinderplus.types._format import (
    AMRFinderPlusAnnotationFormat,
    AMRFinderPlusAnnotationsDirFmt,
    AMRFinderPlusDatabaseDirFmt,
    BinaryFormat,
    TextFormat,
)
from q2_amrfinderplus.types._type import AMRFinderPlusAnnotations, AMRFinderPlusDatabase

citations = Citations.load("citations.bib", package="q2_amrfinderplus")

plugin = Plugin(
    name="amrfinderplus",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-amrfinderplus",
    package="q2_amrfinderplus",
    description="A plugin to find acquired antimicrobial resistance genes and point "
    "mutations in protein and/or assembled nucleotide sequences with "
    "NCBI-AMRFinderPlus.",
    short_description="AMR annotation.",
    citations=[],
)

plugin.methods.register_function(
    function=fetch_amrfinderplus_db,
    inputs={},
    parameters={},
    outputs=[("amrfinderplus_db", AMRFinderPlusDatabase)],
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={
        "amrfinderplus_db": "AMRFinderPlus database.",
    },
    name="Download AMRFinderPlus database.",
    description="Download the latest version of the AMRFinderPlus database.",
    citations=[citations["feldgarden2021amrfinderplus"]],
)

plugin.register_semantic_type_to_format(
    AMRFinderPlusDatabase,
    artifact_format=AMRFinderPlusDatabaseDirFmt,
)
plugin.register_semantic_type_to_format(
    SampleData[AMRFinderPlusAnnotations],
    artifact_format=AMRFinderPlusAnnotationsDirFmt,
)

plugin.register_formats(
    AMRFinderPlusDatabaseDirFmt,
    TextFormat,
    BinaryFormat,
    AMRFinderPlusAnnotationFormat,
    AMRFinderPlusAnnotationsDirFmt,
)
