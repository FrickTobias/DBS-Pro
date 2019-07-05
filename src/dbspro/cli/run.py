"""
Run DPS-Pro pipeline
"""

import logging
import os
import shutil
import sys
import pkg_resources
from snakemake import snakemake

from dbspro.utils import available_cpu_count

logger = logging.getLogger(__name__)


def main(args):
    if not args.targets:
        args.targets = ['umi-counts.txt', 'umi-density-plot.png', 'read-density-plot.png']

    if os.path.isdir(args.directory):
        if args.force and input(f'Directory {args.directory} exists. Remove? (y/n)') == 'y':
            shutil.rmtree(args.directory)
    else:
        os.mkdir(args.directory)

    targets_with_path = [f'{args.directory}/{t}' for t in args.targets]
    # Modified from: https://github.com/NBISweden/IgDiscover/
    snakefile_path = pkg_resources.resource_filename('dbspro', 'Snakefile')
    logger.root.handlers = []
    success = snakemake(snakefile_path,
                        snakemakepath='snakemake',  # Needed in snakemake 3.9.0
                        dryrun=args.dryrun,
                        cores=args.cores,
                        printshellcmds=True,
                        targets=targets_with_path)

    sys.exit(0 if success else 1)


def add_arguments(parser):
    parser.add_argument("-n", "--dryrun", default=False, action='store_true',
                        help="Perform dry run of pipeline.")
    parser.add_argument("-j", "--cores", "--jobs", metavar="<N>", type=int,
                        default=available_cpu_count(),
                        help="Maximum number of cores to run in parallel. Default: Use as manny as available.")
    parser.add_argument('targets', nargs='*', default=[],
                        help='File(s) to create (without paths). If omitted, the full pipeline is run.')
    parser.add_argument('-d', '--directory',
                        help='Path to directory in which to run pipeline and store output. '
                             'Should contain input file (or symbolic link to file).')
    parser.add_argument('-f','--force', default=False, action='store_true',
                        help='Force to run analysis, removes existing analysis if present.')
