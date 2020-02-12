"""
Analyzes demultiplexed and error corrected data and output data files umi_count.tsv and read_count.tsv with
ABC umi and read counts for each DBS. Also output some statistics.
"""

import logging
import pandas as pd
import dnaio
from collections import defaultdict, Counter
from tqdm import tqdm
import os

from dbspro.utils import print_stats

logger = logging.getLogger(__name__)


def main(args):
    # Barcode processing
    logger.info(f"Starting analysis")
    logger.info(f"Saving DBS information to RAM")
    summary = Counter()

    # Set names for ABCs. Creates dict with file names as keys and selected names as values.
    abc_names = {file_name: os.path.basename(file_name).split('-')[0] for file_name in args.umi_abc}

    bc_dict = dict()
    with dnaio.open(args.dbs, mode="r", fileformat="fasta") as reader:
        for read in tqdm(reader):
            bc_dict[read.name] = read.sequence
    logger.info(f"Finished processing DBS sequences")

    # Counting UMI:s found in the different ABC:s for all barcodes.
    logger.info(f"Calculating stats")
    result_dict = dict()
    for current_abc in args.umi_abc:
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
                if bc not in result_dict:
                    result_dict[bc] = {abc_names[abc]: defaultdict(int) for abc in args.umi_abc}

                result_dict[bc][abc_names[current_abc]][read.sequence] += 1

        logger.info(f"Finished reading file: {current_abc}")

    # Barcode-globbing umi/read counter for all ABC:s
    df_umis, df_reads, abc_counter = make_df_from_dict(result_dict, abc_names, sum_filter=args.filter)

    logging.info(f"Writing output files")
    df_umis.to_csv("umi_counts.tsv", sep="\t")
    df_reads.to_csv("read_counts.tsv", sep="\t")

    summary["Total DBS count"] = len(result_dict)

    print_stats(summary, name=__name__)

    # Reporting stats to terminal
    data_to_print = list()
    for abc in abc_names.values():
        data_to_print.append({
            "ABC": abc,
            "Total # UMI": sum(abc_counter[abc]['umis']),
            "N50(UMI/DBS)": n50_counter(abc_counter[abc]['umis']),
            "Total # Reads": sum(abc_counter[abc]['reads']),
            "N50(Reads/DBS)": n50_counter(abc_counter[abc]['reads'])
        })

    print("\nRESULTS")
    df_data = pd.DataFrame(data_to_print).set_index("ABC", drop=True)
    print(df_data)
    print()

    logger.info(f"Finished")


def make_df_from_dict(result_dict, abc_names, sum_filter=0):
    abc_counter = {abc: {"umis": list(), "reads": list()} for abc in abc_names.values()}

    # Output file writing and
    output_umis = list()
    output_reads = list()
    for bc in iter(result_dict):
        output_umis_line = {'BC': bc}
        output_reads_line = output_umis_line.copy()

        for abc in abc_names.values():
            umi_count = len(result_dict[bc][abc])
            read_count = sum(result_dict[bc][abc].values())

            output_umis_line[abc] = umi_count
            output_reads_line[abc] = read_count

            # Add statistics if passing filter.
            if sum(result_dict[bc][abc].values()) >= sum_filter:
                # Add number of UMI:s
                abc_counter[abc]['umis'].append(umi_count)
                # Add number of reads
                abc_counter[abc]['reads'].append(read_count)

        output_umis.append(output_umis_line)
        output_reads.append(output_reads_line)

    # Create dataframe with barcode as index and columns with ABC data.
    df_umis = pd.DataFrame(output_umis, columns=["BC"] + sorted(abc_names.values())).set_index("BC", drop=True)
    df_reads = pd.DataFrame(output_reads, columns=["BC"] + sorted(abc_names.values())).set_index("BC", drop=True)

    return df_umis, df_reads, abc_counter


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


def add_arguments(parser):
    # Arguments
    parser.add_argument("dbs", help="Reads with only DBS seq in fastq format.")
    parser.add_argument("umi_abc", nargs='+',
                        help="Reads with only UMI seq (unique molecular identifier) for ABC (antibody barcode) files"
                             " in fastq format")
    # Options
    parser.add_argument("-f", "--filter", type=int, default=0, help="Number of minimum reads required for an ABC "
                                                                    "to be included in output. DEFAULT: 0")
