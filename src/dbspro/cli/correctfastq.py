"""
Combines starcode output files into raw fastq files creating error corrected fastq files
"""
import logging
from tqdm import tqdm
import dnaio
from collections import Counter
import os
import statistics

from dbspro.utils import print_stats

logger = logging.getLogger(__name__)


def main(args):
    logger.info(f"Starting analysis")
    logger.info(f"Processing file: {args.err_corr}")

    summary = Counter()

    if os.stat(args.err_corr).st_size == 0:
        logging.warning(f"File {args.err_corr} is empty.")

    err_corr = dict()
    reads_per_cluster = list()
    seqs_per_cluster = list()
    for cluster_seq, num_reads, raw_seqs in tqdm(parse_starcode_file(args.err_corr), desc="Parsing clusters"):
        summary["Clusters"] += 1
        reads_per_cluster.append(num_reads)
        seqs_per_cluster.append(len(raw_seqs))
        err_corr.update({raw_seq: cluster_seq for raw_seq in raw_seqs})

    # Add statistics
    summary["Max reads per cluster"] = max(reads_per_cluster)
    summary["Mean reads per cluster"] = statistics.mean(reads_per_cluster)
    summary["Median reads per cluster"] = statistics.median(reads_per_cluster)
    summary["Clusters with one read"] = sum(1 for r in reads_per_cluster if r == 1)

    summary["Max sequences per cluster"] = max(seqs_per_cluster)
    summary["Mean sequences per cluster"] = statistics.mean(seqs_per_cluster)
    summary["Median sequences per cluster"] = statistics.median(seqs_per_cluster)
    summary["Clusters with one sequence"] = sum(1 for r in seqs_per_cluster if r == 1)

    logger.info(f"Correcting sequences and writing to output file.")

    with dnaio.open(args.raw_fastq, mode="r", fileformat="fastq") as reader, \
            dnaio.open(args.corr_fasta, mode="w", fileformat="fasta") as writer:
        for read in tqdm(reader):

            summary["Reads total"] += 1
            if read.sequence in err_corr:
                read.sequence = err_corr[read.sequence]

                writer.write(read)
                summary["Reads corrected"] += 1
            else:
                summary["Reads without corrected sequence"] += 1

    print_stats(summary, name=__name__)

    logger.info(f"Finished")


def parse_starcode_file(filename):
    with open(filename, "r") as file:
        for line in file:
            try:
                cluster_seq, num_reads, raw_seqs_list = line.split()
            except ValueError:
                logging.warning(f"Non-default starcode output line: {line}")
                continue
            raw_seqs = raw_seqs_list.split(",")
            yield cluster_seq, int(num_reads), raw_seqs


def add_arguments(parser):
    parser.add_argument("raw_fastq", help="Fastq file with raw sequences.")
    parser.add_argument("err_corr", help="Starcode default output with error corrected sequences.")
    parser.add_argument("corr_fasta", help="Output file in fasta with error corrected sequences.")
