"""
Correct FASTQ/FASTA with the corrected sequences from starcode clustering
"""
from collections import defaultdict
import logging
import os
import statistics
from pathlib import Path
from typing import Iterator, Tuple, List, Dict

import dnaio
from xopen import xopen

from dbspro.utils import Summary, tqdm

logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument(
        "input", type=Path,
        help="FASTQ/FASTA with uncorrected sequences."
    )
    parser.add_argument(
        "corrections", type=Path,
        help="Starcode output in format, tab-separate entries: <corrected sequnence>, <read count>, <comma-separated"
             "uncorrected sequences>."
    )
    parser.add_argument(
        "-o", "--output-fasta", type=Path,
        help="Output FASTA with corrected sequences."
    )


def main(args):
    run_correctfastq(
        uncorrected_file=args.input,
        corrections_file=args.corrections,
        corrected_fasta=args.output_fasta,
    )


def run_correctfastq(
    uncorrected_file: str,
    corrections_file: str,
    corrected_fasta: str,
):
    logger.info("Starting analysis")
    logger.info(f"Processing file: {corrections_file}")

    summary = Summary()

    if os.stat(corrections_file).st_size == 0:
        logging.warning(f"File {corrections_file} is empty.")

    corr_map = get_corrections(corrections_file, summary)

    logger.info("Correcting sequences and writing to output file.")

    with dnaio.open(uncorrected_file, mode="r") as reader, \
            dnaio.open(corrected_fasta, mode="w") as writer:
        for read in tqdm(reader, desc="Parsing reads"):
            summary["Reads total"] += 1
            if read.sequence in corr_map:
                read.sequence = corr_map[read.sequence]

                writer.write(read)
                summary["Reads corrected"] += 1
            else:
                summary["Reads without corrected sequence"] += 1

    summary.print_stats(name=__name__)

    logger.info("Finished")


def parse_starcode_file(filename: Path) -> Iterator[Tuple[str, int, List[str]]]:
    with xopen(filename, "r") as file:
        for line in file:
            try:
                cluster_seq, num_reads, raw_seqs_list = line.split()
            except ValueError:
                logging.warning(f"Non-default starcode output line: {line}")
                continue
            raw_seqs = raw_seqs_list.split(",")
            yield cluster_seq, int(num_reads), raw_seqs


def get_corrections(corrections_file: Path, summary: Summary) -> Dict[str, str]:
    corr_map = {}
    stats = defaultdict(list)
    for cluster_seq, num_reads, raw_seqs in tqdm(parse_starcode_file(corrections_file), desc="Parsing clusters"):
        summary["Clusters"] += 1
        stats["read"].append(num_reads)
        stats["sequence"].append(len(raw_seqs))
        corr_map.update({raw_seq: cluster_seq for raw_seq in raw_seqs})

    # Add statistics
    for stat, values in stats.items():
        summary[f"Max {stat}s per cluster"] = max(values)
        summary[f"Mean {stat}s per cluster"] = statistics.mean(values)
        summary[f"Median {stat}s per cluster"] = statistics.median(values)
        summary[f"Clusters with one {stat}"] = sum(1 for v in values if v == 1)
    return corr_map
