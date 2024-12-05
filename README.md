# q2-amrfinderplus

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

QIIME 2 plugin to find acquired antimicrobial resistance genes and point mutations in protein and/or assembled nucleotide sequences with AMRFinderPlus.

## Installation
To install _q2-amrfinderplus_, follow the steps described in the [QIIME 2 installation instructions](https://docs.qiime2.org/2024.10/install/native/#qiime-2-pathogenome-distribution) for the pathogenome distribution .


## Functionality
This QIIME 2 plugin contains actions used to annotate protein sequences MAGs and contigs with antimicrobial resistance genes. The underlying tool used is [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/), for details on
the package, please refer to the [AMRFinderPlus documentation](https://github.com/ncbi/amr/wiki)). Checkout the [wiki](https://github.com/bokulich-lab/q2-amrfinderplus/wiki) for usage examples of _q2-amrfinderplus_.

| Action                 | Description                                                                                |
|------------------------|--------------------------------------------------------------------------------------------|
| fetch-amrfinderplus-db | Download AMRFinderPlus database.                                                           |
| annotate               | Annotate protein sequences, MAGs or contigs with antimicrobial resistance gene information.|

## Dev environment
This repository follows the _black_ code style. To make the development slightly easier
there are a couple of pre-commit hooks included here that will ensure that your changes
follow that formatting style. Before you start working on the code, please
install the hooks by executing `make dev` in your conda environment. From then on,
they will be run automatically every time you commit any changes.
