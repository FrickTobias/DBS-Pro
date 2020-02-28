"""
Collect statistics from different log files. Is run in one out of two modes.

1) If run in a DBS-Pro working directory (containing the dbspro.yaml
configfile) single sample mode is run. This only collects data from the current directory.

2) If run one step above a DBS-Pro working directory multi sample mode is run. This looks for possible DBS-Pro
working directories and collect data from all of them, tagging them with the directory name.
"""

from pathlib import Path
from collections import OrderedDict, Counter
import os
import logging
import sys
import contextlib

from dbspro.utils import print_stats

logger = logging.getLogger(__name__)

CONFIG_FILE = "dbspro.yaml"

# TODO add ability to collect data from allready present summary.tsv files in multisample mode.


def main(args):

    logging.info(f"Scanning tree starting from directory {args.directory.absolute()}")

    summary = Summary()

    firstdir = True
    multisample = False
    in_workdir = False

    # Transverse tree in depth-first manner
    for root, dirs, files in os.walk(args.directory):
        if firstdir:
            if CONFIG_FILE not in files:
                logging.info("Running multi-sample mode")
                multisample = True
            else:
                logging.info("Running single-sample mode")
            firstdir = False

        if CONFIG_FILE in files:
            logging.info(f"Looking for log files in directory: {Path(root).name}")
            in_workdir = True

        # Only looks for statistics files in folder "log_files".
        if not root.endswith("log_files"):
            continue

        # Folder must be in tree containing a workdir.
        if not in_workdir:
            continue

        sample = None
        if multisample:
            sample = str(Path(root).parent)

        currentdir = Path(root)

        #
        # Add new file types or modify current ones from here ...
        #

        if "cutadapt-extract-dbs.log" in files:
            parse_cutadapt(currentdir / "cutadapt-extract-dbs.log", summary, sample=sample)

        if "correctfastq-dbs.log" in files:
            parse_dbspro(currentdir / "correctfastq-dbs.log", summary, sample=sample)

        if "cutadapt-extract-abc-umi.log" in files:
            parse_cutadapt(currentdir / "cutadapt-extract-abc-umi.log", summary, sample=sample)

        if "cutadapt-id-abc.log" in files:
            parse_cutadapt(currentdir / "cutadapt-id-abc.log", summary, sample=sample)

        if any(f.startswith("splitcluster") for f in files):
            abc_files = sorted(filter(lambda x: x.startswith("splitcluster"), files))
            for abc_file in abc_files:
                parse_dbspro(currentdir / abc_file, summary, sample=sample)

        if "analyze.log" in files:
            parse_dbspro(currentdir / "analyze.log", summary, sample=sample)

        #
        # ... to here.
        #

        in_workdir = False

    print_stats(summary.counts, name=__name__)

    summary.write(args.output, multisample=multisample)


class Summary:
    """Class to handle all stats"""

    def __init__(self):
        self.data = OrderedDict()
        self.samples = list()
        self.counts = Counter()

    def _add_filetype(self, filetype):
        if filetype not in self.data:
            self.data[filetype] = OrderedDict()
            self.counts["Nr filestypes"] += 1

    def _add_sample(self, sample):
        if sample not in self.samples:
            self.samples.append(sample)
            self.counts["Nr samples"] += 1

    def add_value(self, filetype, parameter, value, sample=None):
        self._add_filetype(filetype)

        if sample:
            self._add_sample(sample)

            if parameter not in self.data[filetype]:
                self.data[filetype][parameter] = dict()
                self.counts["Nr parameters"] += 1

            self.data[filetype][parameter][sample] = value
            self.counts["Nr values"] += 1
        else:
            self.data[filetype][parameter] = value
            self.counts["Nr parameters"] += 1
            self.counts["Nr values"] += 1

    def write(self, output_file, multisample=False, sep=","):
        """Write summary results to output_file"""
        if multisample:
            self._write_multi(output_file, sep=sep)
        else:
            self._write_single(output_file, sep=sep)

    def _write_multi(self, output_file, sep="\t"):
        with open_file(output_file) as output:
            header = ["File", "Parameter"] + self.samples
            print(sep.join(header), file=output)

            for file_name, parameters in self.data.items():
                for parameter, sample_values in parameters.items():
                    line = [file_name, parameter]
                    for sample in self.samples:
                        line.append(sample_values[sample] if sample in sample_values else "NaN")

                    print(sep.join(line), file=output)

    def _write_single(self, output_file, sep="\t"):
        with open_file(output_file) as output:
            header = ["File", "Parameter", "Value"]
            print(sep.join(header), file=output)

            for file_name, stats in self.data.items():
                for parameter, value in stats.items():
                    line = [file_name, parameter, value]
                    print(sep.join(line), file=output)


@contextlib.contextmanager
def open_file(filename=None):
    """Handle writing to file or stdout"""
    # Based on https://stackoverflow.com/questions/17602878/how-to-handle-both-with-open-and-sys-stdout-nicely
    if filename and filename != "-":
        file = open(filename, "w")
    else:
        file = sys.stdout

    try:
        yield file
    finally:
        if file is not sys.stdout:
            file.close()


def parse_cutadapt(path, summary, sample=None):
    """Parse Cutadapt log files. Commonly these will look like below:

    ```
    This is cutadapt 2.8 with Python 3.6.7
    Command line parameters: PARAMETERS
    Processing reads on 4 cores in single-end mode ...
    Finished in 2.19 s (22 us/read; 2.75 M reads/minute).

    === Summary ===

    Total reads processed:                 100,290
    Reads with adapters:                   100,290 (100.0%)
    Reads that were too short:                  95 (0.1%)
    Reads that were too long:                   13 (0.0%)
    Reads written (passing filters):       100,182 (99.9%)

    Total basepairs processed:    14,040,600 bp
    Total written (filtered):      2,003,268 bp (14.3%)

    === Adapter 3 ===
    ...
    ...
    ```
    This script adds the data from the line following `=== Summary ===` until the next line staring with `===` to the
    Summary object.
    """
    logger.info(f"Found Cutadapt log file: {path}")
    with open(path) as file:
        collect = False
        for line in file:
            # Collect stats after first line starting with '===' and stop at next.
            if line.startswith("==="):
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

                # Collect parameter and value
                parameter, value = line.strip().split(":", maxsplit=1)
                value = value.strip().replace(",", "")

                # Collect possible additional with percentage information.
                additional = None
                if " " in value:
                    value, additional = value.split(" ", maxsplit=1)

                summary.add_value(
                    filetype=path.stem,
                    parameter=parameter,
                    value=value,
                    sample=sample
                )

                if additional:
                    summary.add_value(
                        filetype=path.stem,
                        parameter=f"{parameter} (%)",
                        value=additional.strip().split("%")[0].replace("(", ""),
                        sample=sample
                    )


def parse_dbspro(path, summary, sample=None):
    """Parse DBS-Pro core log files. Commonly these will look like below:

    ``
    ...
    ...
    2020-02-26 16:47:22 - splitcluster - INFO: Total time to run: 0.8298838138580322 s
    ===========================================
    STATS SUMMARY - dbspro.cli.splitcluster
    -------------------------------------------
    Reads kept:                             457
    DBS clusters linked to ABC:             356
    Total UMIs:                             370
    Total clustered UMIs:                   366
    ===========================================
    ```

    This script adds the data from the line following `---` until the next line staring with `===` to the
    Summary object.
    """
    logger.info(f"Found DBS-Pro core log file: {path}")
    with open(path) as file:
        collect = False
        for line in file:
            # Collect stats after first line starting with '---' and stop when starts with '==='
            if line.startswith("---"):
                collect = True
                continue
            elif line.startswith("===") and collect:
                break

            if collect:
                # Collect parameter and value
                parameter, value = line.strip().split(":", maxsplit=1)
                value = value.strip().replace(",", "")

                summary.add_value(
                    filetype=path.stem,
                    parameter=parameter,
                    value=value,
                    sample=sample
                )


def add_arguments(parser):
    parser.add_argument("-o", "--output", default="summary_metrics.csv", help="Output CSV file name. Default: %(default)s.")
    parser.add_argument("-d", "--directory", type=Path, default=".",
                        help="Path to directory where to search for log files. Default is current directory "
                             "(%(deafult)s). ")
