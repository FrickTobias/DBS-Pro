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

accepted_file = 'reads.fastq.gz'
accepted_file_ext = '.fastq.gz'

default_handle_file = "construct-info/handles.tsv"
default_abc_file = "construct-info/ABC-sequences.tsv"


def main(args):
    # Check if path to output directory is absolute or make it so.
    if not os.path.isabs(args.directory):
        args.directory = f"{os.getcwd()}/{args.directory}"

    # Check if directory exists. If not make one.
    if not os.path.isdir(args.directory):
        os.mkdir(args.directory)
        logging.info(f'Output directory {args.directory} created.')

    # Append full path to targets.
    targets_with_path = [f"{args.directory}/{t}" for t in args.targets]

    # Check if correct file exists in output directory.
    # If not, try using the input fastq argument if given.
    if args.fastq:
        # Check if path to output directory is absolute or make it so.
        if not os.path.isabs(args.fastq):
            args.fastq = f"{os.getcwd()}/{args.fastq}"

        # Check if file has the correct extension
        if not str(args.fastq).endswith(accepted_file_ext):
            raise FileNotFoundError(f'File {args.fastq} is not accepted input.')

        # Create symbolic link to file in output directory.
        os.symlink(args.fastq, f"{args.directory}/{accepted_file}")
        logging.info('Creating symbolic link for input file in output directory.')

    # Check if tsv files with handle and ABC seqeunces submitted or use default.
    if not args.handles_file:
        dbs_pro_folder = os.path.realpath(__file__).rsplit("/", 4)[0]
        args.handles_file = os.path.join(dbs_pro_folder, default_handle_file)

    if not args.abc_file:
        dbs_pro_folder = os.path.realpath(__file__).rsplit("/", 4)[0]
        args.abc_file = os.path.join(dbs_pro_folder, default_abc_file)

    # Create dict containing the paramaters to be passed to the snakefile.
    configs_dict = {
        'dbs_cluster_dist': args.dbs_cluster_dist,
        'abc_cluster_dist': args.abc_cluster_dist,
        'filter_reads': args.filter_reads,
        'abc_sequences': args.abc_file,
        'handles': args.handles_file
    }

    # Lines below are modified from: https://github.com/NBISweden/IgDiscover/
    snakefile_path = pkg_resources.resource_filename("dbspro", 'Snakefile')
    logger.root.handlers = []
    success = snakemake(snakefile_path,
                        snakemakepath='snakemake',  # Needed in snakemake 3.9.0
                        dryrun=args.dryrun,
                        printdag=args.dag,
                        quiet=False if not args.dag else True,
                        config=configs_dict,
                        cores=args.cores,
                        printshellcmds=True,
                        targets=targets_with_path)

    sys.exit(0 if success else 1)


def add_arguments(parser):
    parser.add_argument("-n", "--dryrun", default=False, action='store_true',
                        help="Perform dry run of pipeline. DEFAULT: False.")
    parser.add_argument("--dag", default=False, action='store_true',
                        help="Print the dag in the graphviz dot language. DEFAULT: False. To det output to pdf file, "
                             "pipe output into dot as follows: '$ dbspro run --dag | dot -Tpdf > dag.pdf'")
    parser.add_argument("-j", "--cores", "--jobs", metavar="<JOBS>", type=int, default=available_cpu_count(),
                        help="Maximum number of cores to run in parallel. DEFAULT: Use as many as available.")
    parser.add_argument('targets', nargs='*', metavar='<TARGETS>',
                        default=['umi-counts.txt', 'umi-density-plot.png', 'read-density-plot.png'],
                        help='File(s) to create excluding paths). If omitted, the full pipeline is run.')
    parser.add_argument('-d', '--directory', default=os.getcwd(), type=str, metavar='<DIRECTORY>',
                        help='Path to directory in which to run pipeline and store output. Should contain input '
                             'fastq file (or symbolic link to file) unless given as argument. DEFAULT: CWD')
    parser.add_argument('-f', '--fastq', default=None, type=str, metavar='<FASTQ>',
                        help="Input fastq file. Should have extension '.fastq.gz'. DEFAULT: None")

    configs = parser.add_argument_group('Pipeline configs')
    configs.add_argument('--handles-file', default=None, type=str, metavar='<TSV>',
                         help="Path to tsv file containing the name and sequence of handles h1, h2 and h3 in the "
                              "construct structure. DEFAULT: Use file handles.tsv in construct-info.")
    configs.add_argument('--abc-file', default=None, type=str, metavar='<TSV>',
                         help="Path to tsv file containing the name and sequence of ABCs. "
                              "DEFAULT: Use file ABC-sequences.tsv in construct-info.")
    configs.add_argument('--dbs-cluster-dist', default=2, type=int, metavar="<DISTANCE>",
                         help="Maximum edit distance to cluster DBS sequences in Starcode. DEFAULT: 2")
    configs.add_argument('--abc-cluster-dist', default=1, type=int, metavar="<DISTANCE>",
                         help="Maximum edit distance to cluster ABC sequences in Starcode. DEFAULT: 1")

    configs.add_argument("--filter-reads", type=int, default=4, metavar="<READS>",
                         help="Minimum reads required for an ABC to be included in analysis output. DEFAULT: 4")
