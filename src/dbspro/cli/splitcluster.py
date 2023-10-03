"""
Split ABC FASTQ with UMIs based on DBS cluster and cluster UMIs for each partion using UMI-tools.
"""
from collections import defaultdict
import logging
from pathlib import Path
from typing import Dict, List, Iterable, Iterator

import dnaio
from dnaio import Sequence
from umi_tools import UMIClusterer

from dbspro.utils import Summary, tqdm

logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument(
        "uncorrected_umi_fastq", type=Path,
        help="Input FASTQ for demultiplexed ABC with uncorrected UMI sequences and corrected DBS sequence in header."
    )
    parser.add_argument(
        "-o", "--output-fasta", default="-", type=Path,
        help="Output FASTA for demultiplexed ABC with corrected UMI sequences. Default: write to stdout"
    )
    parser.add_argument(
        "-t", "--threshold", default=1, type=int,
        help="Edit distance threshold to cluster sequences. Defaulf: %(default)s"
    )
    parser.add_argument(
        "-l", "--length", type=int, required=True,
        help="Required length of UMI sequence. Reads of wrong length are filtered out."
    )
    parser.add_argument(
        "-m", "--method", type=str, default="directional",
        choices=["unique", "percentile", "cluster", "adjacency", "directional"],
        help="Select UMItools clustering method. Defaulf: %(default)s"
    )


def main(args):
    run_splitcluster(
        uncorrected_umis=args.uncorrected_umi_fastq,
        output_fasta=args.output_fasta,
        dist_threshold=args.threshold,
        required_length=args.length,
        clustering_method=args.method,
    )


def run_splitcluster(
    uncorrected_umis: str,
    output_fasta: str,
    dist_threshold: int,
    required_length: int,
    clustering_method: str,
):
    logger.info(f"Filtering reads not of length {required_length} bp.")
    summary = Summary()

    # Read ABC fasta with UMI sequences and save read name and sequence.
    with dnaio.open(uncorrected_umis, mode="r") as file:
        name_to_umi_seq = {read.name: read.sequence for read in file if len(read.sequence) == required_length}

    logger.info(f"Starting clustering of UMIs within each DBS clusters using method: {clustering_method}")
    logger.info(f"Writing corrected reads to {output_fasta}")

    # Set clustering method
    # Based on https://umi-tools.readthedocs.io/en/latest/API.html
    clusterer = UMIClusterer(cluster_method=clustering_method)

    with dnaio.open(uncorrected_umis, mode="r") as reader, \
            dnaio.open(output_fasta, fileformat="fasta", mode="w") as writer:
        
        dbs_umis = defaultdict(list)
        dbs_current = None
        for read in tqdm(reader, desc="Parsing reads"):
            # Get DBS sequence
            dbs = read.name.split(" ")[-1]
            # If new DBS sequence, cluster UMIs and write to output
            if dbs != dbs_current:
                if dbs_current:
                    for read in correct_umis(dbs_umis, clusterer, dist_threshold, summary):
                        writer.write(read)
                dbs_current = dbs
                dbs_umis = defaultdict(list)

            dbs_umis[umi].append(read.sequence)

    summary.print_stats(name=__name__)


AliasType = Dict[str, List[str]]


def correct_umis(umis: AliasType, clusterer: UMIClusterer, threshold: int, summary: Summary) -> Iterator[Sequence]:
    umi_counts = {bytes(umi, encoding='utf-8'): len(reads) for umi, reads in umis.items()}
    summary["Total UMIs"] += len(umi_counts)

    for cluster in clusterer(umi_counts, threshold=threshold):
        summary["Total clustered UMIs"] += 1
        seqs = [seq.decode("utf-8") for seq in cluster]
        canonical_sequnce = seqs[0]

        for seq in seqs:
            for read_name in umis[seq]:
                yield Sequence(read_name, canonical_sequnce)
