"""
Combines starcode output files into raw fastq files creating error corrected fastq files
"""
import logging
import sys
import gzip

from tqdm import tqdm

import DBSpro.utils as DBSpro

logger = logging.getLogger(__name__)

def main(args):

    logger.info("Starting analysis")
    logger.info("Processing file " + args.err_corr)
    generator = DBSpro.FileReader(args.err_corr)
    err_corr = dict()
    for line in tqdm(generator.fileReader()):

        import time
        try: cluster_seq, num_reads, raw_seqs_list = line.split()
        except ValueError:
            print("WARNING! Non-default starcode output line:\t'" + line + "'")
        for raw_seq in raw_seqs_list.split(","):
            if not raw_seq in err_corr:
                err_corr[raw_seq] = cluster_seq

    generator.close()
    logger.info("Error corrected sequenced parsed.")

    logger.info("Correcting sequences and writing to output file.")
    generator = DBSpro.FileReader(args.raw_fastq)
    no_err_corr_seq = int()
    tot_reads = int()
    corr_seqs = int()
    with open(args.corr_fastq, 'w') as openout:
        for read in tqdm(generator.fastqReader()):

            tot_reads += 1
            if read.seq in err_corr:
                cluster_seq = err_corr[read.seq]
                read.seq = cluster_seq
                openout.write(read.fastq_string())
                corr_seqs += 1
            else:
                no_err_corr_seq += 1
    generator.close()

    logger.info("Reads total:\t" + str(tot_reads))
    logger.info("Reads corrected:\t" + str(corr_seqs))
    logger.info("Reads without corrected seq:\t" + str(no_err_corr_seq))
    logger.info("Finished")


def add_arguments(parser):
    parser.add_argument("raw_fastq", help="Fastq file with raw sequences.")
    parser.add_argument("err_corr", help="Starcode default output with error corrected sequences.")
    parser.add_argument("corr_fastq", help="Output file in fastq with error corrected sequences.")

    # Options
    parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                       "Not recommended due to different function "
                                                                       "names in python 2 and 3.")