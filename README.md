[![Build Status](https://travis-ci.org/FrickTobias/DBS-Pro.svg?branch=master)](https://travis-ci.org/FrickTobias/DBS-Pro)

# DBS-Pro Analysis

- [About](#About)
- [Setup](#Setup)
- [Usage](#Usage)
- [Publications](#Publications)

## About

This pipeline analyses data sequencing data from DBS-Pro experiments for protein and PrEST quantification.

## Setup

First, make sure [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) is installed on your system.

1. Clone the git repository.

    ```{bash}
    git clone https://github.com/FrickTobias/DBS-Pro
    ```

2. Move into the git folder and install all dependencies in a conda environment.

    ```{bash}
    cd DBS-Pro
    conda env create -n dbspro -f DBS-Pro/environment.yml
    ```

3. Activate the conda environment.

    ```{bash}
    conda activate dbspro
    ```

4. Install the dbspro package.

    ```{bash}
    pip install .
    ```

For development, please use `pip install -e .[dev]`.

## Usage

Prepare a FASTA with each of the antibody barcodes used in your experiment. The entry name will be used to define the
targets. Also make sure that each sequence is prepended with `^`, this is used for demultiplexing. See the example FASTA below:

```{bash}
>ABC01
^ATGCTG
>ABC02
^GTAGAT
>ABC03
^CTAGCA
```

Use `dbspro init` to create an analysis folder. Provide the FASTA with the antibody barcodes (here named `ABCs.fasta`),
an directory name and one or more FASTQ for the samples.

```{bash}
dbspro init --abc ABCs.fasta <output-folder> <sample1.fastq>
```

If you have several samples you could also provide a CSV file in the line format: `</path/to/sample.fastq>,<sample_name>`.
This enables you to name your samples as you wish. With a CSV the initialization is as follows:

```{bash}
dbspro init --abc ABCs.fasta --sample-csv samples.csv <output-folder>
```

Once the directory has been successfully initialized, moving into the directory

```{bash}
cd <output-folder>
```

and check the current (default) configs using

```{bash}
dbspro config
```

Any changes to the configs should be primaraly be done through the `dbspro config` command to validate the parameters. You can check the construct layout by running `dbspro config --print-construct`. Once the configs are updated you are ready to run the full analysis using this command.

```{bash}
dbspro run
```

For more information on how to run use `dbspro run -h`.

## Publications

Checkout version [v0.1](https://github.com/FrickTobias/DBS-Pro/tree/v0.1) for the pipeline used in:

[Stiller, C., Aghelpasand, H., Frick, T., Westerlund, K., Ahmadian, A., & Eriksson Karlström, A. (2019). Fast and efficient Fc-specific photoaffinity labelling to produce antibody-DNA-conjugates. Bioconjugate chemistry](https://pubs.acs.org/doi/abs/10.1021/acs.bioconjchem.9b00548).
