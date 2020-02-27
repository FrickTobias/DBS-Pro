#! /usr/bin/python3

import sys
import gzip
import logging
import os
import re
import subprocess
import dnaio
import pandas as pd

logger = logging.getLogger(__name__)


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


def available_cpu_count():
    """ Number of available virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling
    userspace-only program
    Taken from: http://stackoverflow.com/a/1006301/715090
    """

    # cpuset
    # cpuset may restrict the number of *available* processors
    try:
        m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$',
                      open('/proc/self/status').read())
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                return res
    except IOError:
        pass

    # Python 2.6+
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        pass

    # https://github.com/giampaolo/psutil
    try:
        import psutil
        return psutil.cpu_count()  # psutil.NUM_CPUS on old versions
    except (ImportError, AttributeError):
        pass

    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

        if res > 0:
            return res
    except (AttributeError, ValueError):
        pass

    # Windows
    try:
        res = int(os.environ['NUMBER_OF_PROCESSORS'])

        if res > 0:
            return res
    except (KeyError, ValueError):
        pass

    # jython
    try:
        from java.lang import Runtime
        runtime = Runtime.getRuntime()
        res = runtime.availableProcessors()
        if res > 0:
            return res
    except ImportError:
        pass

    # BSD
    try:
        sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'],
                                  stdout=subprocess.PIPE)
        scStdout = sysctl.communicate()[0]
        res = int(scStdout)

        if res > 0:
            return res
    except (OSError, ValueError):
        pass

    # Linux
    try:
        res = open('/proc/cpuinfo').read().count('processor\t:')

        if res > 0:
            return res
    except IOError:
        pass

    # Solaris
    try:
        pseudo_devices = os.listdir('/devices/pseudo/')
        res = 0
        for psd in pseudo_devices:
            if re.match(r'^cpuid@[0-9]+$', psd):
                res += 1

        if res > 0:
            return res
    except OSError:
        pass

    # Other UNIXes (heuristic)
    try:
        try:
            dmesg = open('/var/run/dmesg.boot').read()
        except IOError:
            dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
            dmesg = dmesgProcess.communicate()[0]

        res = 0
        while '\ncpu' + str(res) + ':' in dmesg:
            res += 1

        if res > 0:
            return res
    except OSError:
        pass

    raise Exception('Can not determine number of CPUs on this system')


def print_stats(summary, name=None, value_width=15, print_to=sys.stderr):
    """
    Prints stats in nice table with two column for the key and value pairs in summary
    :param summary: collections.Coutner object
    :param name: name of script for header e.g. '__name__'
    :param value_width: width for values column in table
    :param print_to: Where to direct output
    """
    # Get widths for formatting
    max_name_width = max(map(len, summary.keys())) + 1
    width = value_width + max_name_width + 1

    # Header
    print("="*width, file=print_to)
    print(f"STATS SUMMARY - {name}", file=print_to)
    print("-"*width, file=print_to)

    # Print stats in columns
    for name, value in summary.items():
        name += ":"
        if type(value) is int:
            print(f"{name:<{max_name_width}} {value:>{value_width},}", file=print_to)
        elif type(value) is float:
            print(f"{name:<{max_name_width+4}} {value:>{value_width}.3f}", file=print_to)
        else:
            print(f"{name:<{max_name_width}} {value:>{value_width}}", file=print_to)
    print("="*width, file=print_to)
