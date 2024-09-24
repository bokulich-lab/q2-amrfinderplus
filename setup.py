# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

description = "A plugin to find acquired antimicrobial resistance genes and point "
"mutations in protein and/or assembled nucleotide sequences with NCBI-AMRFinderPlus."

setup(
    name="q2-amrfinderplus",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="QIIME 2 development team",
    author_email="rischv@ethz.ch",
    description=description,
    url="https://github.com/bokulich-lab/q2-amrfinderplus",
    entry_points={
        "qiime2.plugins": [
            "q2_amrfinderplus=" "q2_amrfinderplus" ".plugin_setup:plugin"
        ]
    },
    package_data={
        "q2_amrfinderplus": ["citations.bib"],
        "q2_amrfinderplus.tests": ["data/*"],
    },
    zip_safe=False,
)
