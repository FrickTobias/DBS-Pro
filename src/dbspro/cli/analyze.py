"""
Analyzes demultiplexed and error corrected data
"""

import logging
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import dnaio
import copy
from collections import defaultdict
from tqdm import tqdm
import os

logger = logging.getLogger(__name__)


def main(args):
    # Barcode processing
    logger.info(f"Starting analysis")
    logger.info(f"Saving DBS information to RAM")

    # Set names for ABCs. Creates dict with file names as keys and selected names as values.
    default_tsv = os.path.join(os.path.dirname(__file__), "../../../construct-info/ABC-sequences.tsv")
    if args.names == "Use file name":
        abc_names = {file_name: file_name.split('/')[-1] for file_name in args.umi_abc}
    elif args.names:
        abc_names = get_names(args.names, args.umi_abc)
    else:
        abc_names = get_names(default_tsv, args.umi_abc)

    bc_dict = dict()
    with dnaio.open(args.dbs, mode="r", fileformat="fasta") as reader:
        for read in tqdm(reader):
            bc_dict[read.name] = read.sequence
    logger.info(f"Finished processing DBS sequences")

    # Counting UMI:s found in the different ABC:s for all barcodes.
    logger.info(f"Calculating stats")
    result_dict = dict()
    umi_without_proper_bc = int()
    for current_abc in args.umi_abc:
        logger.info(f"Reading file: {current_abc}")

        with dnaio.open(current_abc, mode="r", fileformat="fasta") as reader:
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
                    result_dict[bc] = {abc_names[abc]: defaultdict(int) for abc in args.umi_abc}

                result_dict[bc][abc_names[current_abc]][read.sequence] += 1

        logger.info(f"Finished reading file: {current_abc}")

    # Barcode-globbing umi/read counter for all ABC:s
    df_out, abc_counter = make_df_from_dict(result_dict, abc_names, sum_filter=args.filter)

    logging.info(f"Writing output file to: {args.output}")
    df_out.to_csv(args.output, sep="\t")

    logger.info(f"Total DBS count: {len(result_dict):9,}")
    logger.info(f"UMIs without BC: {umi_without_proper_bc:9,}")

    # Reporting stats to terminal
    data_to_print = list()
    for abc in abc_names.values():
        data_to_print.append({
            "ABC": abc[-25:],
            "Total # UMI": sum(abc_counter[abc]['umis']),
            "N50(UMI/DBS)": n50_counter(abc_counter[abc]['umis']),
            "Total # Reads": sum(abc_counter[abc]['reads']),
            "N50(Reads/DBS)": n50_counter(abc_counter[abc]['reads'])
        })

    print("\nRESULTS")
    df_data = pd.DataFrame(data_to_print).set_index("ABC", drop=True)
    print(df_data)
    print()

    # Plotting
    logger.info(f"Prepping data for plot")

    read_dict_for_plotting, umi_dict_for_plotting = format_data_for_plotting(result_dict)

    plot_density_correlation_matrix(args.read_plot, read_dict_for_plotting, abc_names.values())
    plot_density_correlation_matrix(args.umi_plot, umi_dict_for_plotting, abc_names.values())

    logger.info(f"Finished")


def make_df_from_dict(result_dict, abc_names, sum_filter=0):
    abc_counter = {abc: {"umis": list(), "reads": list()} for abc in abc_names.values()}

    # Output file writing and
    output_list = list()
    for bc in iter(result_dict):
        output_line = {'BC': bc}

        for abc in abc_names.values():
            output_line[abc] = len(result_dict[bc][abc])

            # Add statistics if passing filter.
            if sum(result_dict[bc][abc].values()) >= sum_filter:
                # Add number of UMI:s
                abc_counter[abc]['umis'].append(len(result_dict[bc][abc]))
                # Add number of reads
                abc_counter[abc]['reads'].append(sum(result_dict[bc][abc].values()))

        output_list.append(output_line)

    # Create dataframe with barcode as index and columns with ABC data. Export to tsv.
    df = pd.DataFrame(output_list, columns=["BC"] + sorted(abc_names.values())).set_index("BC", drop=True)
    return df, abc_counter


def get_names(tsv_file, file_names):
    """
    Takes a file tsv_file with columns 'Antibody-target' and 'Barocode-sequence' and list file_names and
    combines them into dict with file_name keys and ABC names as values.
    :param tsv_file: Input .tsv. Must contain columns 'Antibody-target' (with names) and 'Barocode-sequence'
    (with sequences)
    :param file_names: list. File names in list. Names must contain barcodes sequence.
    :return: dict
    """
    df = pd.read_csv(tsv_file, sep="\t")
    df_dict = df.to_dict('records')
    results_dict = dict()
    for record in df_dict:
        # Matches tsv file ABC name with file name by filtering file names for the barcode sequence
        file_name = list(filter(lambda abc: record['Barcode-sequence'] in abc, file_names))[0]
        results_dict[file_name] = record['Antibody-target']
    return results_dict


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
    Plot density correlation matrix.
    :param name:
    :param result_dict:
    :param abc_names:
    :return:
    """
    df_list = [dbs_dict for dbs_dict in result_dict.values()]

    df = pd.DataFrame(df_list, columns=abc_names)
    g = sns.pairplot(df, diag_kind="kde", diag_kws=dict(shade=True, bw=.05, vertical=False))

    for x in range(len(abc_names)):
        for y in range(len(abc_names)):
            g.axes[x, y].set_xlim((0, 50))
            g.axes[x, y].set_ylim((0, 50))

    logger.info(f"Plotting: {name}")

    plt.savefig(name)


def add_arguments(parser):
    # Arguments
    parser.add_argument("dbs", help="Reads with only DBS seq in fastq format.")
    parser.add_argument("output", help="output file")
    parser.add_argument("read_plot", help="Filename for output reads/DBS pair plot (will be .png)")
    parser.add_argument("umi_plot", help="Filename for output UMI:s/DBS pair plot (will be .png)")
    parser.add_argument("umi_abc", nargs='+',
                        help="Reads with only UMI seq (unique molecular identifier) for ABC (antibody barcode) files"
                             " in fastq format")
    # Options
    parser.add_argument("-f", "--filter", type=int, default=0, help="Number of minimum reads required for an ABC "
                                                                    "to be included in output. DEFAULT: 0")
    parser.add_argument("-n", "--names", default="Use file name", nargs="?", metavar="<NAMES>",
                        help="Include tsv file with ABC names and barcodes to use for naming output. If no file is "
                             "given the default file from construct-info/ABC-sequences.tsv is used. If omitted the "
                             "file names are used instead")
