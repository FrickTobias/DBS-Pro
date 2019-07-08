"""
Analyzes demultiplexed and error corrected data
"""

import logging
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import dnaio
import copy
from collections import defaultdict
from tqdm import tqdm

logger = logging.getLogger(__name__)


def main(args):
    # Barcode processing
    logger.info(f"Starting analysis")
    logger.info(f"Saving DBS information to RAM")

    bc_dict = dict()
    with dnaio.open(args.dbs, mode="r", fileformat="fastq") as reader:
        for read in tqdm(reader):
            bc_dict[read.name] = read.sequence
    logger.info(f"Finished processing DBS sequences")

    # Counting UMI:s found in the different ABC:s for all barcodes.
    logger.info(f"Calculating stats")
    result_dict = dict()
    umi_without_proper_bc = int()
    for current_abc in args.umi_abc:
        logger.info(f"Reading file: {current_abc}")

        with dnaio.open(current_abc, mode="r", fileformat="fastq") as reader:
            # Loop over reads in file, where read.seq = umi
            for read in tqdm(reader):

                # Try find UMI
                try:
                    bc = bc_dict[read.name]
                except KeyError:
                    umi_without_proper_bc += 1
                    continue

                # If not dbs in result dict, add it and give it a dictionary for every abc
                if bc not in result_dict:
                    result_dict[bc] = {abc: defaultdict(int) for abc in args.umi_abc}

                result_dict[bc][current_abc][read.sequence] += 1

        logger.info("Finished reading file\t" + str(current_abc))

    # Barcode-globbing umi/read counter for all ABC:s
    abc_counter_umi = {abc: [0] for abc in args.umi_abc}
    abc_counter_read = copy.deepcopy(abc_counter_umi)

    # Output file writing and
    output_list = list()
    for bc in iter(result_dict):
        output_line = {'BC': bc}

        for abc in args.umi_abc:
            output_line[abc] = len(result_dict[bc][abc])

            # Add statistics if passing filter.
            if sum(result_dict[bc][abc].values()) >= args.filter:
                # Add number of UMI:s
                abc_counter_umi[abc].append(len(result_dict[bc][abc]))
                # Add number of reads
                abc_counter_read[abc].append(sum(result_dict[bc][abc].values()))

        output_list.append(output_line)

    df_out = pd.DataFrame(output_list)
    df_out.to_csv(args.output, sep="\t")

    logger.info(f"Tot DBS count: {len(result_dict)}")

    # Reporting stats to terminal
    data_to_print = list()
    for abc in args.umi_abc:
        data_to_print.append({
            "ABC": abc,
            "Total # UMI": sum(abc_counter_umi[abc]),
            "N50(UMI/DBS)": n50_counter(abc_counter_umi[abc]),
            "Total # Reads": sum(abc_counter_read[abc]),
            "N50(Reads/DBS)": n50_counter(abc_counter_read[abc])
        })
    print("\nRESULTS")
    df_data = pd.DataFrame(data_to_print).set_index("ABC", drop=True)
    print(df_data)
    print()

    # Plotting
    logger.info("Prepping data for plot")
    read_dict_for_plotting, umi_dict_for_plotting = format_data_for_plotting(result_dict)
    plot_density_correlation_matrix(args.read_plot, read_dict_for_plotting, args.umi_abc)
    plot_density_correlation_matrix(args.umi_plot, umi_dict_for_plotting, args.umi_abc)
    logger.info("Finished")


def n50_counter(input_list):
    """
    Calculates N50 for a given list
    :param input_list: list with numbers (list)
    :return: N50 (same type as elements of list)
    """
    input_list.sort()
    half_tot = sum(input_list)/2

    current_count = 0
    for num in input_list:
        current_count += num
        if current_count >= half_tot:
            return num


def format_data_for_plotting(result_dict):
    """
    Divides result dict to plot-ready dicts
    :param result_dict: dict[dbs][abc][umi] = read_count
    :return:    dict[dbs][abc] = read_count
                dict[dbs][abc] = umi_count
    """
    read_dict_for_plotting = dict()
    umi_dict_for_plotting = dict()
    for dbs in result_dict.keys():
        read_dict_for_plotting[dbs] = dict()
        umi_dict_for_plotting[dbs] = dict()
        for abc in result_dict[dbs].keys():
            read_dict_for_plotting[dbs][abc] = sum(result_dict[dbs][abc].values())
            umi_dict_for_plotting[dbs][abc] = len(result_dict[dbs][abc].keys())

    return read_dict_for_plotting, umi_dict_for_plotting


def plot_density_correlation_matrix(name, result_dict, abc_names):
    """

    :param x:
    :param y:
    :param z:
    :return:
    """

    df_list = list()
    for dbs_dict in result_dict.values():
        df_list.append(dbs_dict)
    df = pd.DataFrame(df_list, columns=abc_names)
    g = sns.pairplot(df, diag_kind="kde", diag_kws=dict(shade=True, bw=.05, vertical=False))

    for x in range(3):
        for y in range(3):
            g.axes[x,y].set_xlim((0, 50))
            g.axes[x,y].set_ylim((0, 50))
    logger.info("Plotting " + name)

    my_path = os.path.abspath(__file__)
    plt.savefig(name)


def add_arguments(parser):

    # Arguments
    parser.add_argument("dbs", help="Reads with only DBS seq in fastq format.")
    parser.add_argument("output", help="output file")
    parser.add_argument("read_plot", help="Filename for output reads/DBS pair plot (will be .png)")
    parser.add_argument("umi_plot", help="Filename for output UMI:s/DBS pair plot (will be .png)")
    parser.add_argument("umi_abc", nargs='+',
                        help="Reads with only UMI seq (unique molecular identifier) file for ABC (antibody barcode) 1 in fastq format")
    # Options
    parser.add_argument("-f", "--filter", type=int, default=0, help="Number of minimum reads required for an ABC "
                                                                    "to be included in output. DEFAULT: 0")