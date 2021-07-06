[![CI Linux](https://github.com/FrickTobias/DBS-Pro/actions/workflows/ci_linux.yaml/badge.svg?branch=master&event=schedule)](https://github.com/FrickTobias/DBS-Pro/actions/workflows/ci_linux.yaml) [![CI MacOS](https://github.com/FrickTobias/DBS-Pro/actions/workflows/ci_macos.yaml/badge.svg?branch=master&event=schedule)](https://github.com/FrickTobias/DBS-Pro/actions/workflows/ci_macos.yaml) <!-- markdownlint-disable MD041-->

# DBS-Pro Analysis

- [About](#About)
- [Setup](#Setup)
- [Usage](#Usage)
- [Demo](#Demo)
- [Development](#Development)
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
    ```

    For reproducibility the `*.lock` files are used.

    2.1. For OSX use:

    ```{bash}
    conda create --name dbspro --file environment.osx-64.lock
    ```

    2.2. For LINUX use:

    ```{bash}
    conda create --name dbspro --file environment.linux-64.lock
    ```

    2.3. Using flexible dependancies (Not recommended)

    ```{bash}
    conda env create --name dbspro --file environment.yml
    ```

    This option will likely introduce newer versions the softwares
    and depenencies which have not yet been tested.

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

### Standard constructs

The most common construct are included as presets which can be initialized using the `-c/--construct` parameter in `dbspro config`. Currently available constructs include:

#### dbspro_v1

```{bash}
Sequence: 5'-CGATGCTAATCAGATCA BDVHBDVHBDVHBDVHBDVH AAGAGTCAATAGACCATCTAACAGGATTCAGGTA XXXXX NNNNNN TTATATCACGACAAGAG-3'
Name:        |       H1      | |       DBS        | |               H2               | |ABC| |UMI | |       H3      |
Size (bp):   |       17      | |        20        | |               34               | | 5 | | 6  | |       17      |
```

This is the DBS-Pro construct used in the publication [Stiller et al. 2019][1].

#### dbspro_v2

```{bash}
Sequence: 5'-CAGTCTGAGCGGTTCAACAGG BDVHBDVHBDVHBDVHBDVH GCGGTCGTGCTGTATTGTCTCCCACCATGACTAACGCGCTTG XXXXX NNNNNN CACCTGACGCACTGAATACGC-3'
Name:        |         H1        | |       DBS        | |                   H2                   | |ABC| |UMI | |         H3        |
Size (bp):   |         21        | |        20        | |                   42                   | | 5 | | 6  | |         21        |
```

This is the DBS-Pro construct currently in use.

#### pba

```{bash}
Sequence: 5'-NNNNNNNNNNNNNNN ACCTGAGACATCATAATAGCA XXXXX NNNNNN CATTACTAGGAATCACACGCAGAT-3'
Name:        |     DBS     | |         H2        | |ABC| |UMI | |          H3          |
Size (bp):   |      15     | |         21        | | 5 | | 6  | |          24          |
```

This is the construct used in the article [Wu et al. 2019][2] which introduces the Proximity Barcoding Assay (PBA).

## Demo

A short demostration of the pipeline and some downstream analysis is available in the following
[Jupyter Notebook](example/example.ipynb). This can also be used to test that the [conda environment](#Setup) is
properly setup.

## Development

For notes on development see [doc/development](docs/development.rst).

## Publications

Checkout version [v0.1](https://github.com/FrickTobias/DBS-Pro/tree/v0.1) for the pipeline used in:

[Stiller, C., Aghelpasand, H., Frick, T., Westerlund, K., Ahmadian, A., & Eriksson Karlström, A. (2019). Fast and efficient Fc-specific photoaffinity labelling to produce antibody-DNA-conjugates. Bioconjugate chemistry][1].

[1]: https://pubs.acs.org/doi/abs/10.1021/acs.bioconjchem.9b00548 "Stiller et al. 2019"
[2]: https://doi.org/10.1038/s41467-019-11486-1 "Wu et al. 2019"
