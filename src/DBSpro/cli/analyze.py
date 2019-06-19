"""
Analyzes demultiplexed and error corrected data
"""

import logging
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from tqdm import tqdm

import DBSpro.utils as DBSpro


logger = logging.getLogger(__name__)


def main(args):
    # Barcode processing
    logger.info("Starting analysis")
    logger.info("Saving DBS information to RAM")

    generator = DBSpro.FileReader(args.dbs)
    bc_dict = dict()
    for read in tqdm(generator.fastqReader()):
        bc_dict[read.header] = read.seq
    logger.info("Finished processing DBS sequences")

    # Counting UMI:s found in the different ABC:s for all barcodes.
    logger.info("Calculating stats")
    result_dict = dict()
    abc_list = args.umi_abc
    umi_without_proper_bc = int()
    for current_abc in abc_list:
        logger.info("Reading file\t" + str(current_abc))
        generator = DBSpro.FileReader(current_abc)

        # Loop over reads in file, where read.seq = umi
        for read in tqdm(generator.fastqReader()):

            # Try find UMI
            try: bc = bc_dict[read.header]
            except KeyError:
                umi_without_proper_bc += 1
                continue

            # If not dbs in result dict, add it and give it a dictionary for every abc
            if not bc in result_dict:
                result_dict[bc] = dict()
                for abc in abc_list:
                    result_dict[bc][abc] = dict()

            # Add +1 to corresponding UMI sequence
            if not read.seq in result_dict[bc][abc]:
                result_dict[bc][current_abc][read.seq] = int()
            result_dict[bc][current_abc][read.seq] += 1

        logger.info("Finished reading file\t" + str(current_abc))

    # Barcode-globbing umi/read counter for all ABC:s
    abc_counter_umi = dict()
    abc_counter_read = dict()
    for abc in abc_list:
        abc_counter_umi[abc] = [0]
        abc_counter_read[abc] = [0]

    # Output file writing and
    with open(args.output, 'w') as openout:
        for bc in result_dict.keys():
            out_string = str()
            for abc in abc_list:
                # Prepping outstring: umi_count(abc1) + \t + umi_count(abc2) + \t + umi_count(abc3) \n
                out_string += str(len(result_dict[bc][abc])) + '\t'
                # Add number of UMI:s
                if sum(result_dict[bc][abc].values()) >= args.filter:
                    abc_counter_umi[abc].append(len(result_dict[bc][abc].keys()))
                    # Add number of reads
                    abc_counter_read[abc].append(sum(result_dict[bc][abc].values()))
                    # If not enough reads, remove entry
                #else:
                #    del result_dict[bc][abc]

            openout.write(out_string + '\n')

    # Reporting stats to terminal
    print()
    logger.info("Tot DBS count:\t" + str(len(result_dict.keys())))
    for abc in abc_list:
        print()
        logger.info("\t" + abc + " tot umi count:\t" + "{:,}".format(sum(abc_counter_umi[abc])))
        logger.info("\t" + abc + " N50 umi per dbs:\t" + "{:,}".format(n50_counter(abc_counter_umi[abc])))
        logger.info("\t" + abc + " tot read count:\t" + "{:,}".format(sum(abc_counter_read[abc])))
        logger.info("\t" + abc + " N50 read per dbs:\t" + "{:,}".format(n50_counter(abc_counter_read[abc])))
    print()

    # Plotting
    logger.info("Prepping data for plot")
    read_dict_for_plotting, umi_dict_for_plotting = format_data_for_plotting(result_dict)
    plot_density_correlation_matrix(args.read_plot,read_dict_for_plotting, abc_list)
    plot_density_correlation_matrix(args.umi_plot,umi_dict_for_plotting, abc_list)
    logger.info("Finished")

def dict_clearer(dictionary):
    """
    Takes a dictionary and removed any keys which does not have any values
    :param dictionary:
    :return:
    """

    for key, val in dictionary.copy().items():
        if len(val) == 0:
            del dictionary[key]

    return dictionary

def n50_counter(input_list):
    """
    Calculates N50 for a given list
    :param list: list with numbers (list)
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