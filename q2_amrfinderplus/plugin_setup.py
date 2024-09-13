# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Citations, Plugin
from q2_amrfinderplus import __version__

citations = Citations.load("citations.bib", package="q2_amrfinderplus")

plugin = Plugin(
    name="amrfinderplus",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-amrfinderplus",
    package="q2_amrfinderplus",
    description="A plugin to find acquired antimicrobial resistance genes and point mutations in protein and/or assembled nucleotide sequences with NCBI-AMRFinderPlus.",
    short_description="AMR annotation.",
    citations=[]
)
