"""
Combines demultiplexed and error corrected FASTA file and output a aggregated TSV file on the format:

    Barcode Target  UMI ReadCount   Sample

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
    logger.info("Starting analysis")
    summary = Counter()

    # Set names for ABCs. Creates dict with file names as keys and selected names as values.
    abc_names = {file_name: os.path.basename(file_name).split('-')[0].split(".")[-1] for file_name in args.targets}
    sample_name = os.path.basename(args.targets[0]).split(".")[0]

    logger.info("Saving DBS information to RAM")
    bc_dict = get_dbs_headers(args.dbs)
    logger.info("Finished processing DBS sequences")

    # Counting UMI:s found in the different ABC:s for all barcodes.
    logger.info("Calculating stats")
    results = dict()
    for current_abc in args.targets:
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

    df, df_filt = make_dataframes(results, args.filter, sample_name=sample_name)

    output_stats(df_filt, sorted(abc_names.values()))

    logger.info("Sorting data")
    df = df.sort_values(["Barcode", "Target", "UMI"])

    logging.info("Writing output")
    df.to_csv(args.output, sep="\t")

    print_stats(summary, name=__name__)

    logger.info("Finished")


def get_dbs_headers(file):
    dbs = dict()
    with dnaio.open(file, mode="r", fileformat="fasta") as reader:
        for read in tqdm(reader, desc="Parsing DBS reads"):
            dbs[read.name] = read.sequence
    return dbs


def output_stats(df_filt, abcs):
    # Perpare data for analysis
    data_to_print = list()
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


def make_dataframes(results, limit, sample_name):
    # Generate filtered and unfiltered dataframes from data.
    output = list()
    output_filt = list()
    for bc, abcs in tqdm(results.items(), desc="Parsing results"):
        for abc, umis in abcs.items():
            is_ok = True if sum(umis.values()) >= limit else False
            for umi, read_count in umis.items():
                line = {
                    "Barcode": bc,
                    "Target": abc,
                    "UMI": umi,
                    "ReadCount": read_count,
                    "Sample": sample_name
                }
                output.append(line)
                if is_ok:
                    output_filt.append(line)

    # Create dataframe with barcode as index and columns with ABC data.
    cols = ["Barcode", "Target", "UMI", "ReadCount", "Sample"]
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
