"""
Collect statistics from different log files
"""

from pathlib import Path
from collections import OrderedDict
import os
import logging

logger = logging.getLogger(__name__)


def main(args):
    summary = OrderedDict()
    for root, dirs, files in os.walk(args.directory):
        if not root.endswith("log_files"):
            continue

        log_files = list(filter(lambda x: x.endswith(".log"), files))

        currentdir = Path(root)

        if "cutadapt-extract-dbs.log" in log_files:
            parse_cutadapt(currentdir / "cutadapt-extract-dbs.log", summary)

        if "correctfastq-dbs.log" in log_files:
            parse_dbspro(currentdir / "correctfastq-dbs.log", summary)

        if "cutadapt-extract-abc-umi.log" in log_files:
            parse_cutadapt(currentdir / "cutadapt-extract-abc-umi.log", summary)

        if "cutadapt-id-abc.log" in log_files:
            parse_cutadapt(currentdir / "cutadapt-id-abc.log", summary)

        if any(f.startswith("splitcluster") for f in files):
            abc_files = sorted(filter(lambda x: x.startswith("splitcluster"), files))
            for abc_file in abc_files:
                parse_dbspro(currentdir / abc_file, summary)

        if "analyze.log" in log_files:
            parse_dbspro(currentdir / "analyze.log", summary)

    write_summary(summary, args.output)


def write_summary(summary, output_file):
    with open(output_file, "w") as output:
        print(f"File\tParameter\tValue")
        for file_name, stats in summary.items():
            for parameter, value in stats.items():
                print(f"{file_name}\t{parameter}\t{value}", file=output)


def parse_cutadapt(path: Path, summary):
    logger.info(f"Found cutadapt log file: {path.name}")
    stats = OrderedDict()

    with open(path, "r") as file:
        collect = False
        for line in file:
            # Collect stats after first line starting with '=' and stop at next.
            if line.startswith("="):
                if not collect:
                    collect = True
                    continue
                else:
                    break

            if collect:
                # Skip empty lines
                if line.strip() == "":
                    continue

                # Skip lines with basepair level stats.
                if "bp" in line:
                    continue

                parameter, value = line.strip().split(":", maxsplit=1)
                value = value.strip().replace(",", "")

                # Collect possible additional with percentage information.
                additional = None
                if " " in value:
                    value, additional = value.split(" ", maxsplit=1)

                stats[parameter] = int(value)

                if additional:
                    stats[f"{parameter} (%)"] = float(additional.strip().split("%")[0].replace("(", ""))

    summary[path.stem] = stats


def parse_dbspro(path: Path, summary):
    logger.info(f"Found DBS-Pro core log file: {path.name}")
    stats = OrderedDict()

    with open(path, "r") as file:
        collect = False
        for line in file:
            # Collect stats after first line starting with '-' and stop when starts with '='
            if line.startswith("-"):
                collect = True
                continue
            elif line.startswith("=") and collect:
                break

            if collect:
                parameter, value = line.strip().split(":", maxsplit=1)
                value = value.strip().replace(",", "")
                if "." in value:
                    stats[parameter] = float(value)
                else:
                    stats[parameter] = int(value)

    summary[path.stem] = stats


def add_arguments(parser):
    parser.add_argument("-o", "--output", default="summary.tsv", help="Output TSV file name. Default: %(default)s.")
    parser.add_argument("-d", "--directory", type=Path, default=".",
                        help="Path to directory where to search for log files. Default is current directory "
                             "(%(deafult)s). ")
