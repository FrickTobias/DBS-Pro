# iSeq Analysis

This pipeline analyses data sequencing data from iSeq experiments for protein and PrEST quantification.

## Setup

Download the pipeline by first navigating to where you want to put the iSeq Analysis folder and clone the GitHub.

```
git clone https://github.com/FrickTobias/iSeq.git
```
#### Requirements

To run the pipeline make sure you have the following required softwares.

- [cutadapt]()
- [starcode]()
- [pigz]()
- [mail]()

## Useage

To run the pipeline write the following command.

```
bash iSeq_automation.sh -p <processors> -m <email> reads.fq output_folder
```

#### Advanced useage

No advanced useage is currently available. 