[![CI Linux](https://github.com/FrickTobias/DBS-Pro/actions/workflows/ci_linux.yaml/badge.svg?branch=master)](https://github.com/FrickTobias/DBS-Pro/actions/workflows/ci_linux.yaml) [![CI MacOS](https://github.com/FrickTobias/DBS-Pro/actions/workflows/ci_macos.yaml/badge.svg?branch=master)](https://github.com/FrickTobias/DBS-Pro/actions/workflows/ci_macos.yaml) <!-- markdownlint-disable MD041-->

# DBS-Pro Analysis

- [About](#About)
- [Setup](#Setup)
- [Usage](#Usage)
- [Demo](#Demo)
- [Development](#Development)
- [Publications](#Publications)

## About

This pipeline analyses data sequencing data from DBS-Pro experiments for protein and PrEST quantification. The DBS-Pro method uses barcoded antibodies for surface protein quantification in droplets. For example to study [single exosomes][3].

<!-- Image generated using DBS-Pro-testdata-0.4 with command `dbspro run --dag | dot -Tpng -Gdpi=300 > dag.png`.-->
![DBS-Pro pipeline overview](https://user-images.githubusercontent.com/27061883/125053336-47936600-e0a5-11eb-99c4-846bd0f056d7.png)
<p align="center"><i>Overview of DBS-Pro pipeline run on three samples.</i></p>

The pipeline takes input of single end FASTQs with a construct such as those specified in [standard constructs](##Standard-constructs). For each sample the [DBS](#DBS) is extracted (`extract_dbs`) and clustered (`dbs_cluster`) to enable error correction of the DBS sequences (`correct_dbs`). At the same time the [ABC](#ABC) and [UMI](#UMI) are extracted from the same read (`extract_abc_umi`)and then the UMIs are demultiplexed based on their ABC (`demultiplex_abc`). For each ABC the UMIs are grouped by DBS then clustered to correct errors (`umi_cluster`). Finaly the corrected sequences are combined into a read specific DBS, ABC and UMI combination that are tallied to create the final output in the form of a TSV (`integrate`). If there are multiple sampels these are also merged to generate a combined TSV (`merge_data`). A final report is also generated to enable some basic QC of the data. Also see the [demo](/example/example.ipynb) for a step-by-step of a typical workflow.   

<sup><a name="DBS"><b>DBS</b></a>: Droplet Barcode Sequence. Reads sharing this sequence originate from the same droplet.</sup><br/>
<sup><a name="ABC"><b>ABC</b></a>: Antibody Barcodes Sequence. Identifies which antibody was present in the droplet.</sup><br/>
<sup><a name="UMI"><b>UMI</b></a>: Unique Molecular Identifier. Identifies how many antibodies with a particular ABC that was present in the droplet.</sup><br/>

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

This is the DBS-Pro construct used in the publication [Banijamali et al. 2022][3].

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

[Stiller, C., Aghelpasand, H., Frick, T., Westerlund, K., Ahmadian, A., & Eriksson Karlström, A. (2019). *Fast and efficient Fc-specific photoaffinity labelling to produce antibody-DNA-conjugates*. Bioconjugate chemistry][1].

Version [v0.3](https://github.com/FrickTobias/DBS-Pro/tree/v0.3) was used in:

[Banijamali, M., Höjer, P., Nagy, A., Hååg, P., Gomero, E. P., Stiller, C., Kaminskyy, V. O., Ekman, S., Lewensohn, R., Karlström, A. E., Viktorsson, K., & Ahmadian, A. (2022). *Characterizing Single Extracellular Vesicles by Droplet Barcode Sequencing for Protein Analysis*. Journal of Extracellular Vesicles, e12277.][3]

[1]: https://pubs.acs.org/doi/abs/10.1021/acs.bioconjchem.9b00548 "Stiller et al. 2019"
[2]: https://doi.org/10.1038/s41467-019-11486-1 "Wu et al. 2019"
[3]: https://doi.org/10.1002/jev2.12277
