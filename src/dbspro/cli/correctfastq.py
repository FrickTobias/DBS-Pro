"""
Combines starcode output files into raw fastq files creating error corrected fastq files
"""
import logging
from tqdm import tqdm
import dnaio
from collections import Counter
import os

from dbspro.utils import print_stats

logger = logging.getLogger(__name__)


def main(args):
    logger.info(f"Starting analysis")
    logger.info(f"Processing file: {args.err_corr}")

    summary = Counter()

    if os.stat(args.err_corr).st_size == 0:
        logging.warning(f"File {args.err_corr} is empty.")

    err_corr = dict()
    clusters = set()
    with open(args.err_corr, "r") as file:
        for line in tqdm(file):
            try:
                cluster_seq, num_reads, raw_seqs_list = line.split()
            except ValueError:
                logging.warning(f"Non-default starcode output line: {line}")
                continue
            clusters.add(cluster_seq)
            for raw_seq in raw_seqs_list.split(","):
                if raw_seq not in err_corr:
                    err_corr[raw_seq] = cluster_seq

    logger.info(f"Clusters: {len(clusters)}")

    logger.info(f"Error corrected sequenced parsed.")

    logger.info(f"Correcting sequences and writing to output file.")

    with dnaio.open(args.raw_fastq, mode="r", fileformat="fastq") as reader, \
            dnaio.open(args.corr_fasta, mode="w", fileformat="fasta") as openout:
        for read in tqdm(reader):

            summary["Reads total"] += 1
            if read.sequence in err_corr:
                read.sequence = err_corr[read.sequence]

                openout.write(read)
                summary["Reads corrected"] += 1
            else:
                summary["Reads without corrected sequence"] += 1

    print_stats(summary, name=__name__)

    logger.info(f"Finished")


def add_arguments(parser):
    parser.add_argument("raw_fastq", help="Fastq file with raw sequences.")
    parser.add_argument("err_corr", help="Starcode default output with error corrected sequences.")
    parser.add_argument("corr_fasta", help="Output file in fasta with error corrected sequences.")
