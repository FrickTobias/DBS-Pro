"""
Run DPS-Pro pipeline
"""
import logging
import sys
import pkg_resources
from snakemake import snakemake

from dbspro.utils import available_cpu_count

logger = logging.getLogger(__name__)


class SnakemakeError(Exception):
    pass


def add_arguments(parser):
    arg = parser.add_argument
    # Positionals
    arg("targets", nargs="*", metavar="<TARGETS>",
        help="File(s) to create excluding paths). If omitted, the full pipeline is run.")
    # Options
    arg("-n", "--dryrun", default=False, action="store_true",
        help="Perform dry run of pipeline. DEFAULT: False.")
    arg("-j", "--cores", "--jobs", metavar="<JOBS>", type=int, default=available_cpu_count(),
        help="Maximum number of cores to run in parallel. DEFAULT: Use as many as available.")
    arg('--keepgoing', '-k', default=False, action='store_true',
        help='If one job fails, finish the others.')
    arg('--unlock', default=False, action='store_true',
        help='Remove a lock on the working directory.')
    arg("--dir", help="Path to analysis directory. DEFAULT: CWD")

    dags = parser.add_mutually_exclusive_group()
    dags.add_argument(
        "--dag", default=False, action="store_true",
        help="Print the dag in the graphviz dot language. DEFAULT: False. To det output to pdf file, "
             "pipe output into dot as follows: '$ dbspro run --dag | dot -Tpdf > dag.pdf'")
    dags.add_argument(
        "--filegraph", default=False, action='store_true',
        help="Print the file graph showing input/output file from rules in the graphviz dot language (requires "
             "graphviz to be installed). Default: %(default)s. To get output to pdf file, pipe output into dot "
             "as follows: blr run --filegraph | dot -Tpdf > filegraph.pdf")


def main(args):
    targets = args.targets if args.targets else None
    try:
        run(args.dryrun, args.cores, args.keepgoing, args.unlock, args.dag, args.filegraph, targets)
    except SnakemakeError:
        sys.exit(1)
    sys.exit(0)


def run(
        dryrun: bool = False,
        cores: int = 4,
        keepgoing: bool = False,
        unlock: bool = False,
        printdag: bool = False,
        printfilegraph: bool = False,
        targets=None,
        workdir=None,
):
    # snakemake sets up its own logging, and this cannot be easily changed
    # (setting keep_logger=True crashes), so remove our own log handler
    # for now
    logger.root.handlers = []
    snakefile_path = pkg_resources.resource_filename("dbspro", "rules.smk")
    success = snakemake(snakefile_path,
                        snakemakepath="snakemake",  # Needed in snakemake 3.9.0
                        dryrun=dryrun,
                        printdag=printdag,
                        printfilegraph=printfilegraph,
                        quiet=False if not printdag else True,
                        cores=cores,
                        keepgoing=keepgoing,
                        unlock=unlock,
                        printshellcmds=True,
                        targets=targets,
                        workdir=workdir)
    if not success:
        raise SnakemakeError()
