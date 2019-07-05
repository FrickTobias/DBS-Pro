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

    # Modified from: https://github.com/NBISweden/IgDiscover/
    snakefile_path = pkg_resources.resource_filename('dbspro', 'Snakefile')
    logger.root.handlers = []
    success = snakemake(snakefile_path,
                        snakemakepath='snakemake',  # Needed in snakemake 3.9.0
                        dryrun=args.dryrun,
                        cores=args.cores,
                        printshellcmds=True,
                        targets=args.targets if args.targets else None)

    sys.exit(0 if success else 1)


def add_arguments(parser):
    parser.add_argument("-n", "--dryrun", default=False, action='store_true',
                        help="Perform dry run of pipeline.")
    parser.add_argument("-j", "--cores", "--jobs", metavar="<N>", type=int,
                        default=available_cpu_count(),
                        help="Maximum number of cores to run in parallel. Default: Use as manny as available.")
    parser.add_argument('targets', nargs='*', default=[],
                        help='File(s) to create. If omitted, the full pipeline is run.')
