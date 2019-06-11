[![Build Status](https://travis-ci.org/FrickTobias/DBSpro.svg?branch=master)](https://travis-ci.org/FrickTobias/DBSpro)

# DBSpro Analysis

This pipeline analyses data sequencing data from DBSpro experiments for protein and PrEST quantification.

## Setup

Download the pipeline by first navigating to where you want to put the DBSpro Analysis folder and clone the git repository.

```
git clone https://github.com/FrickTobias/DBSpro.git
```
#### Requirements

To run the pipeline make sure you have the following required softwares.

- [cutadapt]()
- [starcode]()
- [pigz]()
- [mail]()

Furthermore some python modules are required to be installed under your python 3 installation.

- [numpy]()
- [pandas]()
- [seaborn]()

## Useage

To run the pipeline write the following command.

```
bash DBSpro_automation.sh -p <processors> -m <email> reads.fq output_folder
```

#### Advanced useage

No advanced useage is currently available. 