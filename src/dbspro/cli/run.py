"""
Run DPS-Pro pipeline
"""

import logging
import os
from pathlib import Path
import sys
import pkg_resources
from snakemake import snakemake

from dbspro.utils import available_cpu_count

logger = logging.getLogger(__name__)

accepted_file = 'reads.fastq.gz'
accepted_file_ext = '.fastq.gz'


def main(args):

    cwd = Path(os.getcwd())

    # Check if path to output directory is absolute or make it so.
    if not os.path.isabs(args.directory):
        args.directory = cwd / args.directory

    # Check if directory exists. If not make one.
    if not os.path.isdir(args.directory):
        os.mkdir(args.directory)
        logging.info(f'Output directory {args.directory} created.')

    # If no targets --> run whole pipeline i.e. use final file names.
    if not args.targets:
        args.targets = ['umi-counts.txt', 'umi-density-plot.png', 'read-density-plot.png']

    # Append full path to targets.
    targets_with_path = [str(args.directory / t) for t in args.targets]

    # Check if correct file exists in output directory.
    # If not, try using the input fastq argument if given.
    if os.path.isfile(args.directory / accepted_file):
        logging.info(f'Found correct file in given directory, continuing analysis.')
    elif args.fastq:
        # Check if path to output directory is absolute or make it so.
        if not os.path.isabs(args.fastq):
            args.fastq = cwd / args.fastq

        # Check if file exists.
        if not os.path.isfile(args.fastq):
            raise FileNotFoundError(f'File {args.fastq} not found.')

        # Check if file has the correct extension
        if not str(args.fastq).endswith(accepted_file_ext):
            raise FileNotFoundError(f'File {args.fastq} is not accepted input.')

        # Create symbolic link to file in output directory.
        os.symlink(args.fastq, args.directory / accepted_file)
        logging.info('Creating symbolic link for input file in output directory.')
    else:
        raise FileNotFoundError(f'Required file not found in folder or given as input.')

    # Lines below are modified from: https://github.com/NBISweden/IgDiscover/
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
    parser.add_argument('targets', nargs='*', default=[], metavar='<TARGETS>',
                        help='File(s) to create excluding paths (If omitted, the full pipeline is run).')
    parser.add_argument('-d', '--directory', default=os.getcwd(), type=Path, metavar='<DIRECTORY>',
                        help='Path to directory in which to run pipeline and store output. Should contain input '
                             'fastq file (or symbolic link to file) unless given as argument. Default: CWD')
    parser.add_argument('-f', '--fastq', default=None, type=str, metavar='<FASTQ>',
                        help="Input fastq file. Should have extension '.fastq.gz'. Default: None")
