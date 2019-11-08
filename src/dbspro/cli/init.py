"""
Create and initialize a new analysis directory.
"""
import logging
import os
import os.path
import sys
from pathlib import Path
from importlib_resources import read_binary

logger = logging.getLogger(__name__)

CONFIGURATION_FILE_NAME = "dbspro.yaml"


def add_arguments(parser):
    parser.add_argument("reads", type=Path, help="Read file (.fastq.gz)")
    parser.add_argument("directory", type=Path, help="New analysis directory to create")


def main(args):
    init(args.directory, args.reads)


def init(directory: Path, reads1: Path):
    if " " in str(directory):
        logger.error("The name of the analysis directory must not contain spaces")
        sys.exit(1)

    create_and_populate_analysis_directory(directory, reads1)

    logger.info(f"Directory {directory} initialized.")
    logger.info(
        'Edit %s/%s, then run "cd %s && blr run" to start the analysis',
        directory,
        CONFIGURATION_FILE_NAME,
        directory,
    )


def create_and_populate_analysis_directory(directory: Path, reads: Path):
    try:
        directory.mkdir()
    except OSError as e:
        logger.error(e)
        sys.exit(1)

    # Write the configuration file
    configuration = read_binary("dbspro", CONFIGURATION_FILE_NAME)
    with (directory / CONFIGURATION_FILE_NAME).open("wb") as f:
        f.write(configuration)

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
    os.symlink(src, os.path.join(dirname, target))
