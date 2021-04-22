"""
Utility functions
"""
from collections import Counter
import logging
import sys

import dnaio
import pandas as pd

logger = logging.getLogger(__name__)

IUPAC_MAP = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'R': {'G', 'A'},
    'Y': {'T', 'C'},
    'M': {'C', 'A'},
    'K': {'G', 'T'},
    'W': {'T', 'A'},
    'S': {'G', 'C'},
    'B': {'G', 'T', 'C'},
    'D': {'G', 'T', 'A'},
    'H': {'C', 'T', 'A'},
    'V': {'G', 'C', 'A'},
    'N': {'G', 'C', 'T', 'A'}}


def get_abcs(abc_fasta_file):
    """
    Helper function to get ABC sequences and names into pandas dataframe
    :param abc_fasta_file:
    :return: dataframe:
    """
    with dnaio.open(abc_fasta_file, fileformat="fasta", mode="r") as abc_fasta:
        abc = pd.DataFrame([{"Sequence": entry.sequence, "Target": entry.name} for entry in abc_fasta])
        abc = abc.set_index("Target", drop=False)

    # Loop over sequences and confirm that they are anchored for cutadapt
    length = None
    for i, row in abc.iterrows():
        assert row['Sequence'].startswith('^'), f"Sequnences in {abc_fasta_file} need to be anchored. " \
                                                f"Add '^' to the start of all ABC sequences."
        if not length:
            length = len(row['Sequence'])
        else:
            assert length == len(row['Sequence']), f"Sequnences in {abc_fasta_file} need of same length. "

    return abc


class Summary(Counter):

    def print_stats(self, name=None, value_width=15, print_to=sys.stderr):
        """
        Prints stats in nice table with two column for the key and value pairs in
        summary
        :param name: name of script for header e.g. '__name__'
        :param value_width: width for values column in table
        :param print_to: Where to direct output. Default: stderr
        """
        # Get widths for formatting
        max_name_width = max(map(len, self.keys()), default=10)
        width = value_width + max_name_width + 1

        # Header
        print("="*width, file=print_to)
        print(f"STATS SUMMARY - {name}", file=print_to)
        print("-"*width, file=print_to)

        # Print stats in columns
        for name, value in self.items():
            value_str = str(value)
        if type(value) is int:
            value_str = f"{value:>{value_width},}"
        elif type(value) is float:
            value_str = f"{value:>{value_width+4},.3f}"

            print(f"{name:<{max_name_width}} {value_str}", file=print_to)
        print("="*width, file=print_to)


def jaccard_index(set1, set2) -> float:
    """Calculate the Jaccard Index metric between two sets"""
    return len(set1 & set2) / len(set1 | set2)
