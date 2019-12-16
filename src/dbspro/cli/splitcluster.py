"""
Split ABC file based on DBS cluster and cluster UMIs for each partion using UMI-tools.
"""
import dnaio
from umi_tools import UMIClusterer
from collections import Counter, defaultdict
import time
import logging

logger = logging.getLogger(__name__)

# Counter object for statistics.
stats = Counter()


def get_umis(file, length=0):
    umis = dict()
    for read in file:
        if length == len(read.sequence):
            umis[read.name] = read.sequence
        else:
            stats['Reads filtered out'] += 1
    return umis


def assign_to_dbs(file, umis):
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
    dbs_umis = defaultdict(dict)
    for read in file:
        # Skip reads that are not to be clustered
        if read.name not in umis:
            continue

        # Get sequences
        dbs = read.sequence
        umi = umis[read.name]

        try:
            dbs_umis[dbs][umi].append(read.name)
        except KeyError:
            dbs_umis[dbs][umi] = [read.name]

        stats["Reads kept"] += 1
    return dbs_umis


def main(args):
    logger.info(f"Filtering reads not of length {args.length} bp.")

    time_start = time.time()

    # Read ABC fasta with UMI sequences and save read name and sequence.
    with dnaio.open(args.abcfile, mode="r") as file:
        umis = get_umis(file, length=args.length)

    time_filtered = time.time()
    logger.info(f"Time for filtering: {time_filtered - time_start} s")
    logger.info(f"Assigning UMIs to DBS clusters")

    with dnaio.open(args.dbsfile, mode="r") as file:
        dbs_umis = assign_to_dbs(file, umis)

    logger.info(f"DBS clusters linked to ABC: {len(dbs_umis)}")

    time_assign = time.time()
    logger.info(f"Time for assigning clusters: {time_assign - time_filtered} s")
    logger.info(f"Starting clustering of UMIs within clusters.")

    # Set clustering method
    # Based on https://umi-tools.readthedocs.io/en/latest/API.html
    clusterer = UMIClusterer(cluster_method='directional')

    with dnaio.open(args.output, fileformat="fasta", mode="w") as output:
        for dbs, umis in dbs_umis.items():
            # Encode each UMI for UMITools and perpare counts
            counts = {bytes(umi, encoding='utf-8'): len(reads) for umi, reads in umis.items()}

            stats["Total UMIs"] += len(counts)

            # Cluster umis
            clustered_umis = clusterer(counts, threshold=args.threshold)

            stats["Total clustered UMIs"] += len(clustered_umis)

            # Loop over clusters and write reads with corrected UMI.
            for cluster in clustered_umis:
                seqs = [seq.decode("utf-8") for seq in cluster]
                canonical_sequnce = seqs[0]

                for seq in seqs:
                    for read_name in umis[seq]:
                        read = dnaio.Sequence(read_name, canonical_sequnce)
                        output.write(read)

    time_end = time.time()
    logger.info(f"Time for clustering: {time_end - time_assign} s")
    logger.info(f"Total time to run: {time_end - time_start} s")

    # Send stats to log
    logger.info(f"Reads filtered out: {stats['Reads filtered out']:,}")
    logger.info(f"Reads kept: {stats['Reads kept']}")
    logger.info(f"Total UMIs: {stats['Total UMIs']}")
    logger.info(f"Total clustered UMIs: {stats['Total clustered UMIs']}")


def add_arguments(parser):
    parser.add_argument("dbsfile",
                        help="Input file with corrected DBS sequences. Could be fasta or fastq")
    parser.add_argument("abcfile",
                        help="Input file with ABC UMI sequences. Could be fasta or fastq")

    parser.add_argument("-o", "--output", default=None, type=str,
                        help="Output cluster results as starcode format as txt file.")
    parser.add_argument("-t", "--threshold", default=1, type=int,
                        help="Edit distance threshold to cluster sequences. DEFAULT=1")
    parser.add_argument("-l", "--length", type=int, required=True,
                        help="Length of barcode sequence. DEFAULT=1")
    parser.add_argument("-m", "--method", type=str, choices=["unique", "percentile", "cluster", "adjacency",
                                                             "directional"],
                        help="Select clustering method. DEFAULT=directional")
