"""
Run DPS-Pro pipeline
"""
import logging
import os
import sys
import pkg_resources
from snakemake import snakemake

from dbspro.utils import available_cpu_count

logger = logging.getLogger(__name__)


def main(args):
    # Create dict containing the paramaters to be passed to the snakefile.
    configs_dict = {
        "dbs_cluster_dist": args.dbs_cluster_dist,
        "abc_cluster_dist": args.abc_cluster_dist,
        "filter_reads": args.filter_reads
    }

    # Lines below are modified from: https://github.com/NBISweden/IgDiscover/
    snakefile_path = pkg_resources.resource_filename("dbspro", "rules.smk")
    logger.root.handlers = []
    success = snakemake(snakefile_path,
                        snakemakepath="snakemake",  # Needed in snakemake 3.9.0
                        dryrun=args.dryrun,
                        printdag=args.dag,
                        quiet=False if not args.dag else True,
                        config=configs_dict,
                        cores=args.cores,
                        printshellcmds=True,
                        targets=args.targets,
                        workdir=args.directory)

    sys.exit(0 if success else 1)


def add_arguments(parser):
    # Positionals
    parser.add_argument("targets", nargs="*", metavar="<TARGETS>",
                        help="File(s) to create excluding paths). If omitted, the full pipeline is run.")
    # Options
    parser.add_argument("-n", "--dryrun", default=False, action="store_true",
                        help="Perform dry run of pipeline. DEFAULT: False.")
    parser.add_argument("--dag", default=False, action="store_true",
                        help="Print the dag in the graphviz dot language. DEFAULT: False. To det output to pdf file, "
                             "pipe output into dot as follows: '$ dbspro run --dag | dot -Tpdf > dag.pdf'")
    parser.add_argument("-j", "--cores", "--jobs", metavar="<JOBS>", type=int, default=available_cpu_count(),
                        help="Maximum number of cores to run in parallel. DEFAULT: Use as many as available.")
    parser.add_argument("-d", "--directory", default=os.getcwd(), type=str, metavar="<DIRECTORY>",
                        help="Path to directory in which to run pipeline and store output. Unless given as argument "
                             "the folder should contain a input fastq file (or symbolic link to file) named "
                             "'reads.fastq.gz'. . DEFAULT: CWD")
    parser.add_argument("-f", "--fastq", default=None, type=str, metavar="<FASTQ>",
                        help="Input fastq file. Should have extension '.fastq.gz'. DEFAULT: None")

    # Configs for snakemake rules
    configs = parser.add_argument_group("Pipeline configs")
    configs.add_argument("--dbs-cluster-dist", default=2, type=int, metavar="<DISTANCE>",
                         help="Maximum edit distance to cluster DBS sequences in Starcode. DEFAULT: 2")
    configs.add_argument("--abc-cluster-dist", default=1, type=int, metavar="<DISTANCE>",
                         help="Maximum edit distance to cluster ABC sequences in Starcode. DEFAULT: 1")

    configs.add_argument("--filter-reads", type=int, default=4, metavar="<READS>",
                         help="Minimum reads required for an ABC to be included in analysis output. DEFAULT: 4")
