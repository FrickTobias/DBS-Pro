"""
Tag FASTQ/FASTA with sequence from matching read by name
"""
from collections import defaultdict
from contextlib import ExitStack
from itertools import islice
import logging
import os
import statistics
from pathlib import Path
from typing import Iterator, Tuple, List, Set, Dict, Optional

import dnaio
from xopen import xopen

from dbspro.utils import Summary, IUPAC_MAP, tqdm

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


def main(args):
    run_correctfastq(
        input=args.input,
        annot=args.annot,
        output=args.output_fasta,
        separator=args.separator,
    )


def run_correctfastq(
    input: str,
    annot: str,
    output: str,
    separator: str,
):
    logger.info("Starting")
    logger.info(f"Processing file: {input}")

    summary = Summary()

    logger.info("Annotating sequences and writing to output file.")
    input_format = determine_filetype(input)
    logger.info(f"Input file format: {input_format}")
    output_format = determine_filetype(output) if str(output) != "-" else input_format
    logger.info(f"Output file format: {output_format}")

    with ExitStack() as stack:
        reader = stack.enter_context(dnaio.open(input, mode="r", fileformat=input_format))
        writer = stack.enter_context(dnaio.open(output, mode="w", fileformat=output_format))
        annotator = stack.enter_context(BufferedFASTAReader(annot))

        for read in tqdm(reader, desc="Parsing reads"):
            summary["Reads total"] += 1
            annot_read = annotator[read.name]
            if annot_read:
                read.name = f"{read.name}{separator}{annot_read.sequence}"
                summary["Reads annotated"] += 1
                writer.write(read)

    summary.print_stats(name=__name__)

    logger.info("Finished")


class BufferedFASTAReader:
    """Read FASTA file and buffer records with same read name"""
    def __init__(self, file):
        self._file = dnaio.open(file, mode="r", fileformat=determine_filetype(file))
        self._iter = iter(self._file)
        self._max_buffer_size = 124
        self._read_buffer = {}

    def __iter__(self):
        return self

    def __getitem__(self, name):
        if name in self._read_buffer:
            return self._read_buffer.pop(name)
        
        for record in islice(self._iter, self._max_buffer_size):
            # If read_name in next pair then parser lines are synced --> drop buffer
            if record.name == name:
                self._read_buffer = {}
                return record
            
            self._read_buffer[record.name] = record

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self._file.close()

def determine_filetype(file):
    # Determine if the file is a FASTQ or FASTA file
    with xopen(file) as f:
        first_line = next(f)
        if first_line.startswith(">"):
            return "fasta"
        elif first_line.startswith("@"):
            return "fastq"
        else:
            raise ValueError(f"File {file} is neither FASTQ or FASTA.")
