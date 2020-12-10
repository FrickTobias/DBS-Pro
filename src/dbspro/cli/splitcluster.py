"""
Split ABC FASTQ with UMIs based on DBS cluster and cluster UMIs for each partion using UMI-tools.
"""
import dnaio
from umi_tools import UMIClusterer
from collections import Counter, defaultdict
import logging
from tqdm import tqdm

from dbspro.utils import print_stats

logger = logging.getLogger(__name__)


def map_names_to_umi_seq(file, summary, length_filter):
    name_to_umi_sequence = dict()
    for read in tqdm(file, desc="Filtering ABC reads"):
        summary["ABC reads total"] += 1
        if length_filter == len(read.sequence):
            name_to_umi_sequence[read.name] = read.sequence
            summary["ABC reads filtered"] += 1

    summary["ABC reads filtered (%)"] = 100 * summary["ABC reads filtered"] / summary["ABC reads total"]
    return name_to_umi_sequence


def assign_to_dbs(file, name_to_umi_seq, summary):
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
    dbs_groups = defaultdict(dict)
    for read in tqdm(file, desc="Assigning DBS"):
        summary["DBS reads"] += 1
        # Skip reads that are not to be clustered
        if read.name not in name_to_umi_seq:
            continue

        # Get sequences
        dbs = read.sequence
        umi = name_to_umi_seq[read.name]

        if umi not in dbs_groups[dbs]:
            dbs_groups[dbs][umi] = list()

        dbs_groups[dbs][umi].append(read.name)

        summary["DBS reads matched to ABCs"] += 1

    if summary["DBS reads"]:
        summary["DBS reads matched to ABCs (%)"] = 100 * summary["DBS reads matched to ABCs"] / summary["DBS reads"]

    return dbs_groups


def main(args):
    logger.info(f"Filtering reads not of length {args.length} bp.")
    summary = Counter()

    # Read ABC fasta with UMI sequences and save read name and sequence.
    with dnaio.open(args.uncorrected_umi_fastq, mode="r") as file:
        name_to_umi_seq = map_names_to_umi_seq(file, summary, length_filter=args.length)

    logger.info("Assigning ABC reads to DBS clusters")

    with dnaio.open(args.corrected_dbs_fasta, mode="r") as file:
        dbs_groups = assign_to_dbs(file, name_to_umi_seq, summary)

    summary["DBS clusters linked to ABC"] = len(dbs_groups)

    logger.info(f"Starting clustering of UMIs within each DBS clusters using method: {args.method}")
    logger.info(f"Writing corrected reads to {args.output_fasta}")

    # Set clustering method
    # Based on https://umi-tools.readthedocs.io/en/latest/API.html
    clusterer = UMIClusterer(cluster_method=args.method)

    with dnaio.open(args.output_fasta, fileformat="fasta", mode="w") as output:
        for dbs, umis in tqdm(dbs_groups.items(), desc="Correcting UMIs"):
            # Encode each UMI for UMITools and perpare counts
            counts = {bytes(umi, encoding='utf-8'): len(reads) for umi, reads in umis.items()}

            summary["Total UMIs"] += len(counts)

            # Cluster umis
            clustered_umis = clusterer(counts, threshold=args.threshold)

            summary["Total clustered UMIs"] += len(clustered_umis)

            # Loop over clusters and write reads with corrected UMI.
            for cluster in clustered_umis:
                seqs = [seq.decode("utf-8") for seq in cluster]
                canonical_sequnce = seqs[0]

                for seq in seqs:
                    for read_name in umis[seq]:
                        read = dnaio.Sequence(read_name, canonical_sequnce)
                        output.write(read)

    print_stats(summary, name=__name__)


def add_arguments(parser):
    parser.add_argument("corrected_dbs_fasta",
                        help="Input FASTA with corrected DBS sequences")
    parser.add_argument("uncorrected_umi_fastq",
                        help="Input FASTQ for demultiplexed ABC with uncorrected UMI sequences.")

    parser.add_argument("-o", "--output-fasta", default="-", type=str,
                        help="Output FASTA for demultiplexed ABC with corrected UMI sequences. Default: write to "
                             "stdout")
    parser.add_argument("-t", "--threshold", default=1, type=int,
                        help="Edit distance threshold to cluster sequences. Defaulf: %(default)s")
    parser.add_argument("-l", "--length", type=int, required=True,
                        help="Required length of UMI sequence. Reads of wrong length are filtered out.")
    parser.add_argument("-m", "--method", type=str, default="directional",
                        choices=["unique", "percentile", "cluster", "adjacency", "directional"],
                        help="Select UMItools clustering method. Defaulf: %(default)s")
