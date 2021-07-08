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
        "dbs_file", type=Path,
        help="Path to FASTA with error-corrected DBS barcode sequences."
    )
    parser.add_argument(
        "target_files", nargs="+", type=Path,
        help="Path to ABC-specific FASTAs with UMI sequences to combine with DBS. "
    )
    parser.add_argument(
        "-o", "--output", default=sys.stdout, type=Path,
        help="Output TSV to file instead of stdout."
    )


def main(args):
    run_analysis(
        dbs_file=args.dbs_file,
        target_files=args.target_files,
        output=args.output,
    )


def run_analysis(
    dbs_file: str,
    target_files: List[str],
    output: str,
):
    # Barcode processing
    logger.info("Starting analysis")
    summary = Summary()

    # Set names for ABCs. Creates dict with file names as keys and selected names as values.
    target_file_to_name = {file: os.path.basename(file).split('-')[0].split(".")[-1] for file in target_files}
    sample_name = os.path.basename(target_files[0]).split(".")[0]
    logging.info(f"Found sample {sample_name}.")

    logger.info("Saving DBS information to RAM")
    header_to_dbs = map_header_to_sequence(dbs_file)
    summary["Total DBS reads"] = len(header_to_dbs)

    # Counting UMI:s found in the different ABC:s for all barcodes.
    logger.info("Calculating stats")
    results = get_results(target_files, target_file_to_name, header_to_dbs, summary)

    summary["Total DBS count"] = len(results)

    df = make_dataframe(results)

    # Attach sample info
    df["Sample"] = sample_name

    logger.info("Sorting data")
    df = df.sort_values(["Barcode", "Target", "UMI"])

    logging.info("Writing output")
    df.to_csv(output, sep="\t")

    summary.print_stats(name=__name__)

    logger.info("Finished")


def map_header_to_sequence(file: Path) -> Dict[str, str]:
    with dnaio.open(file, mode="r", fileformat="fasta") as reader:
        return {r.name: r.sequence for r in tqdm(reader, desc="Parsing DBS reads")}


def get_results(target_files: List[Path], target_file_to_name: Dict[Path, str], header_to_dbs: Dict[str, str],
                summary: Dict[str, int]) -> Dict[Tuple[str, str, str], int]:
    results = defaultdict(int)
    for current_target in target_files:
        logger.info(f"Reading file: {current_target}")

        with dnaio.open(current_target, mode="r", fileformat="fasta") as reader:
            # Loop over reads in file, where read.seq = umi
            for read in tqdm(reader, desc=f"Parsing {target_file_to_name[current_target]} reads"):
                summary["Total target reads"] += 1
                dbs = header_to_dbs.get(read.name)

                if dbs is None:
                    summary["Target reads without DBS"] += 1
                    continue

                results[(dbs, target_file_to_name[current_target], read.sequence)] += 1

        logger.info(f"Finished reading file: {current_target}")

    return results


def make_dataframe(results: Dict[Tuple[str, str, str], int]) -> pd.DataFrame:
    output = [(*dbs_target_umi, count) for dbs_target_umi, count in tqdm(results.items(), desc="Parsing results")]
    # Create dataframe with barcode as index and columns with ABC data.
    cols = ["Barcode", "Target", "UMI", "ReadCount"]
    return pd.DataFrame(output, columns=cols).set_index("Barcode", drop=True)
