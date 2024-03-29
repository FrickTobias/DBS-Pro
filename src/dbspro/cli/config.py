"""
Update configuration file. If no --set option is given the current settings are printed.
"""
# Copied & modified from BLR github 16/12 - 2019, https://github.com/FrickTobias/BLR
# Script is based on repos NBISSweden/IgDisover config script.
# Link https://github.com/NBISweden/IgDiscover/blob/master/src/igdiscover/cli/config.py
import sys
import os
import logging
from importlib.resources import files, as_file
from typing import List, Tuple, Dict
from pathlib import Path

from ruamel.yaml import YAML
from snakemake.utils import validate

from dbspro.utils import get_abcs

logger = logging.getLogger(__name__)
DEFAULT_PATH = Path("dbspro.yaml")
SCHEMA_FILE = "config.schema.yaml"

STANDARD_CONSTRUCTS = {
    "dbspro_v1": [
        ("h1", "CGATGCTAATCAGATCA"),
        ("h2", "AAGAGTCAATAGACCATCTAACAGGATTCAGGTA"),
        ("h3", "TTATATCACGACAAGAG"),
        ("dbs", "BDVHBDVHBDVHBDVHBDVH"),
    ],
    "dbspro_v2": [
        ("h1", "CAGTCTGAGCGGTTCAACAGG"),
        ("h2", "GCGGTCGTGCTGTATTGTCTCCCACCATGACTAACGCGCTTG"),
        ("h3", "CACCTGACGCACTGAATACGC"),
        ("dbs", "BDVHBDVHBDVHBDVHBDVH"),
    ],
    "pba": [
        ("h1", "null"),
        ("h2", "ACCTGAGACATCATAATAGCA"),
        ("h3", "CATTACTAGGAATCACACGCAGAT"),
        ("dbs", "NNNNNNNNNNNNNNN"),
    ],
    "dbspro_v3": [
        ("h1", "GGCGTTACAGCATGGATGTGG"),
        ("h2", "GCGCAGTTACATGAGACTCTGCCTGGCCGCGATTGCAGATGTTA"),
        ("h3", "TGGCGTGGAGTGCAGGTATAC"),
        ("dbs", "BDVHBDVHBDVHBDVHBDVH"),
    ],
}


def main(args):
    run_config(args.file, args.set, args.construct, args.print_construct)


def run_config(file: Path, change_set: List[Tuple] = None, construct: str = None, print_construct: bool = False):
    if change_set is None:
        change_set = []

    if change_set or construct is not None:
        if construct is not None:
            change_set.extend(STANDARD_CONSTRUCTS[construct])
        change_config(file, change_set)
    elif print_construct:
        print_construct(file)
    else:
        print_config(file)


def print_config(file: Path):
    configs, yaml = load_yaml(file)
    print(f"--- CONFIGS IN: {file} ---")
    yaml.dump(configs, stream=sys.stdout)


def print_construct(file: Path):
    configs, yaml = load_yaml(file)

    # Get handles
    h1 = Handle("H1", configs["h1"]) if configs["h1"] is not None else Handle("H1", "")
    dbs = Handle("DBS", configs["dbs"])
    h2 = Handle("H2", configs["h2"])
    umi = Handle("UMI", "N"*configs["umi_len"])
    abcs = get_abcs(configs["abc_file"])
    abc = Handle("ABC", "X"*len(abcs["Sequence"].iloc[0].strip("^")))
    h3 = Handle("H3", configs["h3"])

    handles = [h1, dbs, h2, abc, umi, h3] if configs["h1"] is not None else [dbs, h2, abc, umi, h3]

    print(
        f"--- CONSTRUCT LAYOUT (Total length = {sum(h.length for h in handles)}) ---\n"
        f"Sequence: 5'-{' '.join(h.seq for h in handles)}-3'\n"
        f"Name:        {' '.join(h.name for h in handles)}\n"
        f"Size (bp):   {' '.join(h.size for h in handles)}\n"
    )


class Handle:
    def __init__(self, name, sequence):
        self.seq = sequence
        self.length = len(self.seq)
        self.size = self._add_padding(str(self.length))
        self.name = self._add_padding(name)

    def _add_padding(self, string: str):
        return f"|{string.center(self.length-2, ' ')}|"


def change_config(filename: Path, changes_set: List[Tuple[str, str]]):
    """
    Change config YAML file at filename using the changes_set key-value pairs.
    :param filename: string with path to YAML config file to change.
    :param changes_set: dict with changes to incorporate.
    """
    # Get configs from file.
    configs, yaml = load_yaml(filename)

    # Update configs
    for key, value in changes_set:
        update_configs(configs, key, value)

    # Confirm that configs is valid.
    with as_file(files("dbspro").joinpath(SCHEMA_FILE)) as schema_path:
        validate(configs, schema_path)

    # Write first to temporary file then overwrite filename.
    tmpfile = str(filename) + ".tmp"
    with open(tmpfile, "w") as file:
        yaml.dump(configs, stream=file)
    os.rename(tmpfile, filename)


def update_configs(configs: Dict[str, str], key: str, value: str):
    """
    Check key and value before updating configs
    """
    key = key.lower()

    if key not in configs:
        logger.warning(f"KEY = {key} not in config. Config not updated with set ({key}, {value})")
        return

    value = YAML(typ='safe').load(value)
    logger.info(f"Changing value of '{key}': {configs[key]} --> {value}.")
    configs[key] = value


def load_yaml(filename: Path) -> Tuple[Dict[str, str], YAML]:
    """
    Load YAML file and return the yaml object and data.
    :param filename: Path to YAML file
    :return: (data, yaml).
    """
    with open(filename) as file:
        yaml = YAML()
        data = yaml.load(file)
    return data, yaml


def add_arguments(parser):
    parser.add_argument(
        "-s", "--set", nargs=2, metavar=("KEY", "VALUE"), action="append", default=[],
        help="Set KEY to VALUE. Use KEY.SUBKEY[.SUBSUBKEY...] for nested keys. For empty values "
             "write 'null'. Can be given multiple times."
    )
    parser.add_argument(
        "--file", default=DEFAULT_PATH, type=Path,
        help="Configuration file to modify. Default: %(default)s in current directory."
    )
    parser.add_argument(
        "-p", "--print-construct", default=False, action="store_true",
        help="Show construct layout with current handles inplace."
    )
    parser.add_argument(
        "-c", "--construct", choices=list(STANDARD_CONSTRUCTS),
        help="Setup one of the standard contructs."
    )
