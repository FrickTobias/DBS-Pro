[![Build Status](https://travis-ci.org/FrickTobias/DBS-Pro.svg?branch=master)](https://travis-ci.org/FrickTobias/DBS-Pro)

# DBS-Pro Analysis

This pipeline analyses data sequencing data from DBS-Pro experiments for protein and PrEST quantification.

## Setup

First, make sure [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) is installed on your system.

1. Clone the git repository.

```
git clone https://github.com/FrickTobias/DBS-Pro
```

2. Move into the git folder and install all dependencies in a conda environment.

```
cd DBS-Pro
conda env create -n dbspro -f DBS-Pro/environment.yml
```

3. Activate the conda environment.

```
conda activate dbspro
```

4. Install the dbspro package.

```
pip install -e .
```

## Usage

Create an analysis folder with the required input files.
```
dbspro init <raw-reads.fastq> <output-folder>
```
Run the analysis.
```
cd <output-folder>
dbspro run -t <threads>
```
For more information use the `-h` (`--help`) option. 

## Publications

[Stiller, C., Aghelpasand, H., Frick, T., Westerlund, K., Ahmadian, A., & Eriksson Karlström, A. (2019). Fast and efficient Fc-specific photoaffinity labelling to produce antibody-DNA-conjugates. Bioconjugate chemistry.](https://pubs.acs.org/doi/abs/10.1021/acs.bioconjchem.9b00548).
