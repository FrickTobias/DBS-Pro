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
from typing import List, Iterator, Tuple

logger = logging.getLogger(__name__)

CONFIGURATION_FILE_NAME = "dbspro.yaml"
ABC_FILE_NAME = "ABC-sequences.fasta"
ACCEPTED_FILE_EXT = ".fastq.gz"
SAMPLE_FILE_NAME = "samples.tsv"


def add_arguments(parser):
    parser.add_argument(
        "directory", type=Path, help="New analysis directory to create"
    )
    parser.add_argument(
        "reads", nargs="*", type=Path,
        help="Read file(s) (.fastq.gz). Can also be supplied through '-s/--sample-csv'"
    )
    parser.add_argument(
        "-s", "--sample-csv", type=Path,
        help="Path to CSV with one sample per line. Line fromat: <path/to/sample.fastq.gz>,<sample_name>."
    )
    parser.add_argument(
        "--abc", required=True, type=Path, metavar="ABC-sequences.fasta",
        help="Antibody barcode (ABC) sequence fasta file. Should contain the target name in the "
             "header and the ABC seqeunce for demuliplexing."
    )
    return parser


def main(args):
    init(args.directory, args.reads, args.abc, args.sample_csv)


def init(directory: Path, reads: List[Path], abc: Path, sample_csv: str = None):
    if " " in str(directory):
        logger.error("The name of the analysis directory must not contain spaces")
        sys.exit(1)

    if reads != [] and sample_csv is not None:
        logger.error("Only provide a sample CSV or paths to samples, not both.")
        sys.exit(1)

    if reads == [] and sample_csv is None:
        logger.error("Provide a sample CSV or paths to reads.")
        sys.exit(1)

    create_and_populate_analysis_directory(directory, reads, abc, sample_csv)

    logger.info(f"Directory {directory} initialized.")
    logger.info(
        "Edit %s/%s and run 'cd %s && dbspro run' to start the analysis.",
        directory,
        CONFIGURATION_FILE_NAME,
        directory
    )


def create_and_populate_analysis_directory(directory: Path, reads: List[Path], abc_file: Path, sample_csv: Path):
    try:
        directory.mkdir()
    except OSError as e:
        logger.error(e)
        sys.exit(1)

    # Write the configuration file
    configuration = read_binary("dbspro", CONFIGURATION_FILE_NAME)
    with (directory / CONFIGURATION_FILE_NAME).open("wb") as f:
        f.write(configuration)

    # Write ABC FASTA after checking that its correctly formatted
    write_abc_fasta_to_dir(abc_file, directory)

    # Symlink sample FASTQs into workdir and create TSV with sample info
    with (directory / SAMPLE_FILE_NAME).open(mode="w") as f:
        print("Sample", "Reads", "FastqPath", sep="\t", file=f)
        for file, name in get_path_and_name(reads, sample_csv):
            logger.info(f"File {file.name} given sample name '{name}'")
            create_symlink(file, directory, name + ".fastq.gz")
            print(name, count_reads(file), file.resolve(), sep="\t", file=f)


def get_path_and_name(reads: List[Path], sample_csv: Path) -> Iterator[Tuple[Path, str]]:
    if reads != []:
        for file in reads:
            name = fix_name(file.name)
            yield file, name
    else:
        with sample_csv.open() as f:
            for line in f:
                file, name = line.strip().split(",", maxsplit=1)
                name = fix_name(name)
                file = Path(file)
                assert file.exists()
                yield file, name


def fix_name(name: str) -> str:
    return name.replace(".fastq.gz", "").replace("-", "_").replace(".", "_")


def fail_if_inaccessible(path: Path):
    try:
        with path.open():
            pass
    except OSError as e:
        logger.error("Could not open %r: %s", path, e)
        sys.exit(1)


def create_symlink(readspath: Path, dirname: Path, target: str):
    if not os.path.isabs(readspath):
        src = os.path.relpath(readspath, dirname)
    else:
        src = readspath

    # Check if file has the correct extension
    if not str(readspath).endswith(ACCEPTED_FILE_EXT):
        raise FileNotFoundError(f"File {readspath} is not accepted input.")

    os.symlink(src, os.path.join(dirname, target))


def count_reads(fastq: Path) -> int:
    with dnaio.open(fastq, mode="r") as f:
        return sum([1 for _ in f])


def write_abc_fasta_to_dir(abc_file: Path, directory: Path):
    abc_file_out = str(directory / ABC_FILE_NAME)
    with dnaio.open(abc_file_out, mode='w', fileformat="fasta") as open_out, \
            dnaio.open(abc_file, mode='r', fileformat="fasta") as open_in:
        for read in open_in:
            if "." in read.name or "-" in read.name:
                new_name = read.name.replace(".", "_").replace("-", "_")
                logging.info(f"Renaming target {read.name} -> {new_name}")
                read.name = new_name

            if not read.sequence.startswith("^"):
                read.sequence = "^" + read.sequence

            open_out.write(read)
