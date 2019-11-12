"""
Create and initialize a new analysis directory.
"""
import logging
import os
import os.path
import sys
import dnaio
from pathlib import Path
from importlib_resources import read_binary
import pkg_resources

logger = logging.getLogger(__name__)

CONFIGURATION_FILE_NAME = "dbspro.yaml"
ABC_FILE_NAME = "ABC-sequences.fasta"
ACCEPTED_FILE_EXT = ".fastq.gz"


def add_arguments(parser):
    parser.add_argument("reads", type=Path, help="Read file (.fastq.gz)")
    parser.add_argument("directory", type=Path, help="New analysis directory to create")
    parser.add_argument("--abc", default=pkg_resources.resource_filename("dbspro", "ABC-sequences.fasta"), type=str,
                        metavar="ABC-sequences.fasta",
                        help="Antibody barcode (ABC) sequence fasta file")

    return parser


def main(args):
    init(args.directory, args.reads, args.abc)



def init(directory: Path, reads: Path, abc: str):
    if " " in str(directory):
        logger.error("The name of the analysis directory must not contain spaces")
        sys.exit(1)

    create_and_populate_analysis_directory(directory, reads, abc)

    logger.info(f"Directory {directory} initialized.")
    logger.info(
        "Edit %s/%s and run 'cd %s && dbspro run' to start the analysis.",
        directory,
        CONFIGURATION_FILE_NAME,
        directory
    )


def create_and_populate_analysis_directory(directory: Path, reads: Path, abc: str):
    try:
        directory.mkdir()
    except OSError as e:
        logger.error(e)
        sys.exit(1)

    # Write the configuration file
    configuration = read_binary("dbspro", CONFIGURATION_FILE_NAME)
    with (directory / CONFIGURATION_FILE_NAME).open("wb") as f:
        f.write(configuration)

    # Write ABC fasta file with ^ prior to sequence (used in cutadapt command)
    filename_as_string = str(directory) + "/" + ABC_FILE_NAME
    with dnaio.open(filename_as_string, mode='w', fileformat="fasta") as open_out, \
            dnaio.open(abc, mode='r', fileformat="fasta") as open_in:
        for abc in open_in:
            if not abc.sequence.startswith("^"):
                abc.sequence = "^" + abc.sequence

            open_out.write(abc)

    create_symlink(reads, directory, "reads.fastq.gz")


def fail_if_inaccessible(path):
    try:
        with path.open():
            pass
    except OSError as e:
        logger.error("Could not open %r: %s", path, e)
        sys.exit(1)


def create_symlink(readspath, dirname, target):
    if not os.path.isabs(readspath):
        src = os.path.relpath(readspath, dirname)
    else:
        src = readspath

    # Check if file has the correct extension
    if not str(readspath).endswith(ACCEPTED_FILE_EXT):
        raise FileNotFoundError(f"File {readspath} is not accepted input.")

    os.symlink(src, os.path.join(dirname, target))
