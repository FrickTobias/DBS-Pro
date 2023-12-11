"""
Tag FASTQ/FASTA with sequence from matching read by name
"""
from contextlib import ExitStack
from itertools import islice
import logging

from pathlib import Path

import dnaio

from dbspro.utils import Summary, tqdm

logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument(
        "input", type=Path,
        help="Input FASTQ/FASTA to annotate."
    )
    parser.add_argument(
        "annot", type=Path,
        help="FASTQ/FASTA used to annotate input."
    )
    parser.add_argument(
        "-o", "--output-fasta", type=Path, default="-",
        help="Output FASTA with corrected sequences."
    )
    parser.add_argument(
        "-s", "--separator", default="_",
        help="Separetor used to connect annotation string to read name."
    )
    parser.add_argument(
        "-b", "--buffer-size", type=int, default=64,
        help="Buffer size for annotation file."
    )


def main(args):
    run_tagfastq(
        input=args.input,
        annot=args.annot,
        output=args.output_fasta,
        separator=args.separator,
        buffer_size=args.buffer_size,
    )


def run_tagfastq(
    input: str,
    annot: str,
    output: str,
    separator: str,
    buffer_size: int,
):
    logger.info("Starting")
    logger.info(f"Processing file: {input}")
    logger.info(f"Annotating with file: {annot}")

    summary = Summary()

    logger.info("Annotating sequences and writing to output file.")
    input_format = determine_filetype(input)
    logger.info(f"Input file format: {input_format}")
    output_format = determine_filetype(output) if str(output) not in {"-", "/dev/null"} else input_format
    logger.info(f"Output file format: {output_format}")

    with ExitStack() as stack:
        reader = stack.enter_context(dnaio.open(input, mode="r", fileformat=input_format))
        writer = stack.enter_context(dnaio.open(output, mode="w", fileformat=output_format))
        annotator = stack.enter_context(BufferedFASTAReader(annot, buffer_size=buffer_size))

        for read in tqdm(reader, desc="Parsing reads"):
            summary["Reads total"] += 1
            sequence = annotator.get(read.id)
            if sequence:
                read.name = f"{read.id}{separator}{sequence}"
                summary["Reads annotated"] += 1
                writer.write(read)

    summary.print_stats(name=__name__)

    logger.info("Finished")


class BufferedFASTAReader:
    """Read FASTA file and buffer records with same read name"""
    __slots__ = ["_file", "_iter", "_names", "_seqs", "_buffer_size", "_missed", "_count"]

    def __init__(self, file, buffer_size: int = 64):
        self._file = dnaio.open(file, mode="r", fileformat=determine_filetype(file))
        self._buffer_size = buffer_size
        self._iter = iter(self._file)
        self._names = []
        self._seqs = []
        self._missed = 0
        self._count = 0

    def __iter__(self):
        return self

    def get(self, name):
        try:
            index = self._names.index(name)
        except ValueError:
            for record in islice(self._iter, self._buffer_size-len(self._names)):
                self._count += 1
                # If read_name in next pair then parser lines are synced --> drop buffer)
                if record.id == name:
                    self._names.clear()
                    self._seqs.clear()
                    self._missed = 0
                    return record.sequence

                self._names.append(record.id)
                self._seqs.append(record.sequence)

            self._missed += 1

        else:
            sequence = self._seqs[index]
            self._names = self._names[index+1:]
            self._seqs = self._seqs[index+1:]
            self._missed = 0
            return sequence

        if self._missed > self._buffer_size:
            self._buffer_size *= 2
            logger.warning(f"Resetting buffer after {self._missed} misses buffer = {self._buffer_size}.")
            self._missed = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self._file.close()


def determine_filetype(file):
    # Determine if the file is a FASTQ or FASTA file
    file_name = str(file).strip(".gz")
    if file_name.endswith(".fasta") or file_name.endswith(".fa"):
        return "fasta"
    elif file_name.endswith(".fastq") or file_name.endswith(".fq"):
        return "fastq"
    raise ValueError(f"File {file} is neither FASTQ or FASTA.")
