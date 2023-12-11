"""
Run DPS-Pro pipeline
"""
import logging
import subprocess
import sys
from importlib.resources import files, as_file
from typing import List, Optional
from pathlib import Path

from snakemake.utils import available_cpu_count

logger = logging.getLogger(__name__)


class SnakemakeError(Exception):
    pass


def add_arguments(parser):
    arg = parser.add_argument
    # Options
    arg("-c", "--cores", type=int, default=available_cpu_count(),
        help="Number of cores to use for Snakemake. Default: %(default)s (all available cores).")
    arg("--no-use-conda", default=False, action="store_true",
        help="Skip passing argument '--use-conda' to snakemake.")

    # This argument will not capture any arguments due to nargs=-1. Instead parse_known_args()
    # is used in __main__.py to add any arguments not captured here to snakemake_args.
    smk_args = parser.add_argument_group("snakemake arguments")
    smk_args.add_argument(
        'snakemake_args', nargs=-1,
        help="Arguments passed to snakemake. For info about snakemake options run "
             "'snakemake --help'."
    )


def main(args):
    try:
        run(
            cores=args.cores,
            no_conda=args.no_use_conda,
            snakemake_args=args.snakemake_args
        )
    except SnakemakeError:
        sys.exit(1)
    sys.exit(0)


def run(
        cores: int = 1,
        no_conda: bool = False,
        workdir: Optional[Path] = None,
        snakemake_args: Optional[List[str]] = None,
):
    # snakemake sets up its own logging, and this cannot be easily changed
    # (setting keep_logger=True crashes), so remove our own log handler
    # for now
    logger.root.handlers = []
    with as_file(files("dbspro").joinpath("rules.smk")) as snakefile_path:
        cmd = ["snakemake", "-s", str(snakefile_path)]
        cmd += ["--cores", str(cores)]

        # Set defaults
        cmd += ["--printshellcmds"]

        if not no_conda:
            cmd += ["--use-conda"]

        if workdir is not None:
            cmd += ["--directory", str(workdir)]

        if snakemake_args is not None:
            cmd += snakemake_args

        logger.debug(f"Command: {' '.join(cmd)}")
        subprocess.check_call(cmd)
