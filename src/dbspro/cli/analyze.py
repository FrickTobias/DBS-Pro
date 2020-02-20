"""
Combines demultiplexed and error corrected FASTA file and output a aggregated TSV file on the format:

    Barcode TargetName  UMI ReadCount

Output statistics based on filter.

"""

import logging
import pandas as pd
import dnaio
from collections import defaultdict, Counter
from tqdm import tqdm
import os
import sys

from dbspro.utils import print_stats

logger = logging.getLogger(__name__)


def main(args):
    # Barcode processing
    logger.info(f"Starting analysis")
    summary = Counter()

    # Set names for ABCs. Creates dict with file names as keys and selected names as values.
    abc_names = {file_name: os.path.basename(file_name).split('-')[0] for file_name in args.targets}

    logger.info(f"Saving DBS information to RAM")
    bc_dict = get_dbs_headers(args.dbs)
    logger.info(f"Finished processing DBS sequences")

    # Counting UMI:s found in the different ABC:s for all barcodes.
    logger.info(f"Calculating stats")
    results = dict()
    for current_abc in args.targets:
        logger.info(f"Reading file: {current_abc}")

        with dnaio.open(current_abc, mode="r", fileformat="fasta") as reader:
            # Loop over reads in file, where read.seq = umi
            for read in tqdm(reader):

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

    df, df_filt = make_dataframes(results, args.filter)

    output_stats(df_filt, sorted(abc_names.values()))

    logger.info("Sorting data")
    df = df.sort_values(["Barcode", "TargetName", "UMI"])

    logging.info(f"Writing output")
    df.to_csv(args.output, sep="\t")

    print_stats(summary, name=__name__)

    logger.info(f"Finished")


def get_dbs_headers(file):
    dbs = dict()
    with dnaio.open(file, mode="r", fileformat="fasta") as reader:
        for read in tqdm(reader):
            dbs[read.name] = read.sequence
    return dbs


def output_stats(df_filt, abcs):
    # Perpare data for analysis
    data_to_print = list()
    for abc in abcs:
        logger.info(f"Adding data for {abc}")
        data = df_filt[df_filt.TargetName.eq(abc)].groupby("Barcode")
        umis = data.count()["UMI"].tolist()
        reads = data.sum()["ReadCount"].tolist()

        data_to_print.append({
            "ABC": abc,
            "Total # UMI": sum(umis),
            "N50(UMI/DBS)": n50_counter(umis),
            "Total # Reads": sum(reads),
            "N50(Reads/DBS)": n50_counter(reads)
        })

    print("\nRESULTS")
    cols = ["ABC", "Total # UMI", "N50(UMI/DBS)", "Total # Reads", "N50(Reads/DBS)"]
    df_data = pd.DataFrame(data_to_print, columns=cols).set_index("ABC", drop=True)
    print(df_data)
    print()


def n50_counter(input_list):
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


def make_dataframes(results, limit):
    # Generate filtered and unfiltered dataframes from data.
    output = list()
    output_filt = list()
    for bc, abcs in tqdm(results.items(), desc="Parsing results"):
        for abc, umis in abcs.items():
            is_ok = True if sum(umis.values()) >= limit else False
            for umi, read_count in umis.items():
                line = {
                    "Barcode": bc,
                    "TargetName": abc,
                    "UMI": umi,
                    "ReadCount": read_count
                }
                output.append(line)
                if is_ok:
                    output_filt.append(line)

    # Create dataframe with barcode as index and columns with ABC data.
    cols = ["Barcode", "TargetName", "UMI", "ReadCount"]
    return pd.DataFrame(output, columns=cols).set_index("Barcode", drop=True), \
           pd.DataFrame(output_filt, columns=cols).set_index("Barcode", drop=True)


def add_arguments(parser):
    # Arguments
    parser.add_argument("dbs",
                        help="Path to FASTA with error-corrected DBS barcode sequences.")
    parser.add_argument("targets", nargs="+",
                        help="Path to directory containing error-corrected target (ABC) FASTAs with UMI "
                             "(unique molecular identifier) sequences.")

    parser.add_argument("-o", "--output", default=sys.stdout, help="Output TSV. Defualt: %(default)s")

    parser.add_argument("-f", "--filter", type=int, default=0, help="Number of minimum reads required for an ABC "
                                                                    "to be included in output. DEFAULT: 0")
