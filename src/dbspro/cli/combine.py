"""
Combine sequences of two FASTQ/FASTA files so that output is CSV with 'file1_seq, file2_seq' pairs.
"""
import sys
import dnaio
import contextlib


def main(args):
    with dnaio.open(args.file1, mode="r") as reads:
        header_sequences = parse_and_build_dict(reads)

    with dnaio.open(args.file2, mode="r") as reads, open_out(args.output) as out:

        for read in reads:
            try:
                paired_seq = header_sequences[read.name]
            except KeyError:
                continue

            print(f"{paired_seq},{read.sequence}", file=out)


def parse_and_build_dict(file):
    header_sequences = dict()
    for read in file:
        header_sequences[read.name] = read.sequence
    return header_sequences


@contextlib.contextmanager
def open_out(file_name):
    if file_name is sys.stdout:
        yield file_name
    else:
        yield open(file_name, "w")


def add_arguments(parser):
    parser.add_argument("file1", help="FASTQ/FASTA file 1")
    parser.add_argument("file2", help="FASTQ/FASTA file 2")
    parser.add_argument("-o", "--output", help="Output file name. Default: stdout", default=sys.stdout)
