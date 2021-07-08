"""
Split ABC FASTQ with UMIs based on DBS cluster and cluster UMIs for each partion using UMI-tools.
"""
from collections import defaultdict
import logging
from pathlib import Path
from typing import Dict, List, Iterable, Iterator

import dnaio
from dnaio import Sequence
from tqdm import tqdm
from umi_tools import UMIClusterer

from dbspro.utils import Summary

logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument(
        "corrected_dbs_fasta", type=Path,
        help="Input FASTA with corrected DBS sequences"
    )
    parser.add_argument(
        "uncorrected_umi_fastq", type=Path,
        help="Input FASTQ for demultiplexed ABC with uncorrected UMI sequences."
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
        corrected_barcodes=args.corrected_dbs_fasta,
        uncorrected_umis=args.uncorrected_umi_fastq,
        output_fasta=args.output_fasta,
        dist_threshold=args.threshold,
        required_length=args.length,
        clustering_method=args.method,
    )


def run_splitcluster(
    corrected_barcodes: str,
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

    logger.info("Assigning ABC reads to DBS clusters")

    with dnaio.open(corrected_barcodes, mode="r") as file:
        dbs_groups = assign_to_dbs(file, name_to_umi_seq, summary)

    summary["DBS clusters linked to ABC"] = len(dbs_groups)

    logger.info(f"Starting clustering of UMIs within each DBS clusters using method: {clustering_method}")
    logger.info(f"Writing corrected reads to {output_fasta}")

    # Set clustering method
    # Based on https://umi-tools.readthedocs.io/en/latest/API.html
    clusterer = UMIClusterer(cluster_method=clustering_method)

    with dnaio.open(output_fasta, fileformat="fasta", mode="w") as output:
        for _, umis in tqdm(dbs_groups.items(), desc="Correcting UMIs"):
            for read in correct_umis(umis, clusterer, dist_threshold, summary):
                output.write(read)

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


def assign_to_dbs(file: Iterable[Sequence], name_to_umi: Dict[str, str], summary: Summary) -> Dict[str, AliasType]:
    """
    Create structure for storing DBS, UMI and read name

           {DBS1:{UMI1:[name1, name2, ...],
                  UMI2:[name3, ...],
                  ...
                  }
            DBS2:{UMI1:[name1, name2, ...],
                  UMI2:[name3, ...],
                  ...
                  }
           }
    """
    dbs_groups = defaultdict(lambda: defaultdict(list))
    for read in tqdm(file, desc="Assigning DBS"):
        summary["DBS reads"] += 1
        # Get sequences
        dbs = read.sequence
        umi = name_to_umi.get(read.name)

        if umi is None:
            continue

        dbs_groups[dbs][umi].append(read.name)

        summary["DBS reads matched to ABCs"] += 1

    if summary["DBS reads"]:
        summary["DBS reads matched to ABCs (%)"] = 100 * summary["DBS reads matched to ABCs"] / summary["DBS reads"]

    return dbs_groups
