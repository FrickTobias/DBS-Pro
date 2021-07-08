"""
Combines DBS and ABC-demultiplexed UMI FASTA file(s) into TSV file.

Each TSV row has the following format:

    Barcode Target  UMI ReadCount   Sample
"""

import logging
from collections import defaultdict
import os
import sys
from typing import List, Dict, Tuple
from pathlib import Path

import dnaio
import pandas as pd
from tqdm import tqdm

from dbspro.utils import Summary

logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument(
        "dbs", type=Path,
        help="Path to FASTA with error-corrected DBS barcode sequences."
    )
    parser.add_argument(
        "targets", nargs="+", type=Path,
        help="Path to ABC-specific FASTAs with UMI sequences to combine with DBS. "
    )
    parser.add_argument(
        "-o", "--output", default=sys.stdout, type=Path,
        help="Output TSV to file instead of stdout."
    )


def main(args):
    run_analysis(
        dbs=args.dbs,
        targets=args.targets,
        output=args.output,
    )


def run_analysis(
    dbs: str,
    targets: List[str],
    output: str,
):
    # Barcode processing
    logger.info("Starting analysis")
    summary = Summary()

    # Set names for ABCs. Creates dict with file names as keys and selected names as values.
    abc_names = {file_name: os.path.basename(file_name).split('-')[0].split(".")[-1] for file_name in targets}
    sample_name = os.path.basename(targets[0]).split(".")[0]

    logger.info("Saving DBS information to RAM")
    bc_dict = get_dbs_headers(dbs)
    logger.info("Finished processing DBS sequences")

    # Counting UMI:s found in the different ABC:s for all barcodes.
    logger.info("Calculating stats")
    results = dict()
    for current_abc in targets:
        logger.info(f"Reading file: {current_abc}")

        with dnaio.open(current_abc, mode="r", fileformat="fasta") as reader:
            # Loop over reads in file, where read.seq = umi
            for read in tqdm(reader, desc=f"Parsing {abc_names[current_abc]} reads"):

                # Try find UMI
                try:
                    bc = bc_dict[read.name]
                except KeyError:
                    summary["UMIs without BC"] += 1
                    continue

                # If not dbs in result dict, add it and give it a dictionary for every abc
                if bc not in results:
                    results[bc] = {abc: defaultdict(int) for abc in abc_names.values()}

                results[bc][abc_names[current_abc]][read.sequence] += 1

        logger.info(f"Finished reading file: {current_abc}")

    summary["Total DBS count"] = len(results)

    df = make_dataframe(results, sample_name=sample_name)

    logger.info("Sorting data")
    df = df.sort_values(["Barcode", "Target", "UMI"])

    logging.info("Writing output")
    df.to_csv(output, sep="\t")

    summary.print_stats(name=__name__)

    logger.info("Finished")


def get_dbs_headers(file: Path) -> Dict[str, str]:
    with dnaio.open(file, mode="r", fileformat="fasta") as reader:
        return {r.name: r.sequence for r in tqdm(reader, desc="Parsing DBS reads")}


AliasType = Dict[str, Dict[str, Dict[str, Dict[str, int]]]]


def make_dataframe(results: AliasType, sample_name: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # Generate filtered and unfiltered dataframes from data.
    output = list()
    for bc, abcs in tqdm(results.items(), desc="Parsing results"):
        for abc, umis in abcs.items():
            for umi, read_count in umis.items():
                line = {
                    "Barcode": bc,
                    "Target": abc,
                    "UMI": umi,
                    "ReadCount": read_count,
                    "Sample": sample_name
                }
                output.append(line)

    # Create dataframe with barcode as index and columns with ABC data.
    cols = ["Barcode", "Target", "UMI", "ReadCount", "Sample"]
    return pd.DataFrame(output, columns=cols).set_index("Barcode", drop=True)
