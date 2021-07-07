"""
Combines demultiplexed and error corrected FASTA file and output a aggregated TSV file on the format:

    Barcode Target  UMI ReadCount   Sample

Output statistics based on filter.

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
        help="Path to directory containing error-corrected target (ABC) FASTAs with UMI "
             "(unique molecular identifier) sequences."
    )
    parser.add_argument(
        "-o", "--output", default=sys.stdout, type=Path,
        help="Output TSV. Defualt: %(default)s"
    )
    parser.add_argument(
        "-f", "--filter", type=int, default=0,
        help="Number of minimum reads required for an ABC to be included in output. Default: %(default)s"
    )


def main(args):
    run_analysis(
        dbs_file=args.dbs,
        target_files=args.targets,
        output=args.output,
        filter=args.filter,
    )


def run_analysis(
    dbs_file: str,
    target_files: List[str],
    output: str,
    filter: int
):
    # Barcode processing
    logger.info("Starting analysis")
    summary = Summary()

    # Set names for ABCs. Creates dict with file names as keys and selected names as values.
    target_file_to_name = {file: os.path.basename(file).split('-')[0].split(".")[-1] for file in target_files}
    sample_name = os.path.basename(target_files[0]).split(".")[0]

    logger.info("Saving DBS information to RAM")
    header_to_dbs = map_header_to_sequence(dbs_file)
    logger.info("Finished processing DBS sequences")

    # Counting UMI:s found in the different ABC:s for all barcodes.
    logger.info("Calculating stats")
    results = get_results(target_files, target_file_to_name, header_to_dbs, summary)

    summary["Total DBS count"] = len(results)

    df, df_filt = make_dataframes(results, filter, sample_name=sample_name)

    output_stats(df_filt, sorted(target_file_to_name.values()))

    logger.info("Sorting data")
    df = df.sort_values(["Barcode", "Target", "UMI"])

    logging.info("Writing output")
    df.to_csv(output, sep="\t")

    summary.print_stats(name=__name__)

    logger.info("Finished")


def get_results(target_files, file_to_name, header_to_dbs, summary):
    results = {}
    for target_file in target_files:
        target_name = file_to_name[target_file]
        logger.info(f"Reading file: {target_file}")

        with dnaio.open(target_file, mode="r", fileformat="fasta") as reader:
            # Loop over reads in file, where read.seq = umi
            for read in tqdm(reader, desc=f"Parsing {file_to_name[target_file]} reads"):

                dbs = header_to_dbs.get(read.name)

                if dbs is None:
                    summary["UMIs without BC"] += 1

                # If not dbs in result dict, add it and give it a dictionary for every abc
                if dbs not in results:
                    results[dbs] = {abc: defaultdict(int) for abc in file_to_name.values()}

                results[dbs][target_name][read.sequence] += 1

        logger.info(f"Finished reading file: {target_file}")
    return results


def map_header_to_sequence(file: Path) -> Dict[str, str]:
    with dnaio.open(file, mode="r", fileformat="fasta") as reader:
        return {r.name: r.sequence for r in tqdm(reader, desc="Parsing DBS reads")}


def output_stats(df_filt: pd.DataFrame, abcs: List[str]):
    # Perpare data for analysis
    data_to_print = []
    for abc in abcs:
        logger.info(f"Adding data for {abc}")
        data = df_filt[df_filt.Target.eq(abc)].groupby("Barcode")
        umis = data.count()["UMI"].tolist()
        reads = data.sum()["ReadCount"].tolist()

        data_to_print.append({
            "ABC": abc,
            "Total # UMI": sum(umis),
            "N50(UMI/DBS)": n50_counter(umis),
            "Total # Reads": sum(reads),
            "N50(Reads/DBS)": n50_counter(reads)
        })

    sys.stderr.write("\nRESULTS\n")
    cols = ["ABC", "Total # UMI", "N50(UMI/DBS)", "Total # Reads", "N50(Reads/DBS)"]
    df_data = pd.DataFrame(data_to_print, columns=cols).set_index("ABC", drop=True)
    sys.stderr.write(df_data.to_string())
    sys.stderr.write("\n")


def n50_counter(input_list: List[int]) -> int:
    """
    Calculates N50 for a given list
    :param input_list: list with numbers (list)
    :return: N50 (same type as elements of list)
    """
    input_list.sort()
    half_tot = sum(input_list) / 2

    current_count = 0
    for num in input_list:
        current_count += num
        if current_count >= half_tot:
            return num


AliasType = Dict[str, Dict[str, Dict[str, Dict[str, int]]]]


def make_dataframes(results: AliasType, count_threshold: int, sample_name: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # Generate filtered and unfiltered dataframes from data.
    output = []
    output_filt = []
    for dbs, abc_to_umis in tqdm(results.items(), desc="Parsing results"):
        for abc, umi_to_count in abc_to_umis.items():
            passed_threshold = sum(umi_to_count.values()) >= count_threshold

            for umi, read_count in umi_to_count.items():
                line = {
                    "Barcode": dbs,
                    "Target": abc,
                    "UMI": umi,
                    "ReadCount": read_count,
                    "Sample": sample_name
                }
                output.append(line)
                if passed_threshold:
                    output_filt.append(line)

    # Create dataframe with barcode as index and columns with ABC data.
    cols = ["Barcode", "Target", "UMI", "ReadCount", "Sample"]
    output = pd.DataFrame(output, columns=cols).set_index("Barcode", drop=True)
    output_filt = pd.DataFrame(output_filt, columns=cols).set_index("Barcode", drop=True)

    return output, output_filt
