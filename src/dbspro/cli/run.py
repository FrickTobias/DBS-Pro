"""
Run DPS-Pro pipeline
"""
import logging
import sys
import pkg_resources
from snakemake import snakemake

from dbspro.utils import available_cpu_count

logger = logging.getLogger(__name__)


def main(args):
    # Lines below are modified from: https://github.com/NBISweden/IgDiscover/
    snakefile_path = pkg_resources.resource_filename("dbspro", "rules.smk")
    logger.root.handlers = []
    success = snakemake(snakefile_path,
                        snakemakepath="snakemake",  # Needed in snakemake 3.9.0
                        dryrun=args.dryrun,
                        printdag=args.dag,
                        quiet=False if not args.dag else True,
                        cores=args.cores,
                        printshellcmds=True,
                        targets=args.targets,
                        workdir=args.dir)

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
    parser.add_argument("--dir", help="Path to analysis directory. DEFAULT: CWD")
