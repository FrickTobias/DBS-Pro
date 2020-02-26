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


class FileReader(object):
    """
    Reads input files as generator, handles gzip.
    """

    def __init__(self, filehandle, filehandle2=None):

        """
        Setup function, detects if files are gzipped and saves file handles (generator =
        FileReader(filehandle=args.input_file)). If only one file is to be read, only use first the first argument.
        :param filehandle: string. File handle name. Typically args.input_file.
        :param filehandle2: string OR None. Second file handle name, if only one file should be read leave blank.
        """
        # Init variables setting
        self.filehandle = filehandle
        self.gzip = bool()

        if self.filehandle == "stdin":
            self.openfile = sys.stdin
        # Open files as zipped or not not (depending on if they end with .gz)
        elif self.filehandle[-3:] == '.gz':
            logger.info('File detected as gzipped, unzipping when reading')
            self.openfile = gzip.open(self.filehandle, 'r')
            self.gzip = True
        else:
            self.openfile = open(self.filehandle, 'r')

        # Paired end preparation
        self.filehandle2 = filehandle2
        if self.filehandle2:

            # Open files as zipped or not not (depending on if they end with .gz)
            if self.filehandle2[-3:] == '.gz':
                logger.info('File detected as gzipped, unzipping when reading')

                self.openfile2 = gzip.open(self.filehandle2, 'r')
            else:
                self.openfile2 = open(self.filehandle2, 'r')

    def fileReader(self):
        """
        Reads non-specific (non-structured) files as generator.
        :return: strin. Yields one line for every iteration.
        """
        for line in self.openfile:
            if self.gzip:
                line = line.decode("utf-8")
            yield line

    def fastqReader(self):
        """
        Reads fastq format files as generator, reads 4 lines at the time (=one read).
        :return: instance. Fastq reads as instances (see BLR FastqRead object function).
        """

        line_chunk = list()
        for line in self.openfile:
            if self.gzip:
                line = line.decode("utf-8")
            line_chunk.append(line)
            if len(line_chunk) == 4:
                read = FastqRead(line_chunk)
                line_chunk = list()
                yield read

    def fastqPairedReader(self):
        """
        Reads two paired fastq files as generator and yields a pair of two reads.
        :return: instance, instance. read1 and read2 as instances (see BLR FastqRead object function).
        """

        line_chunk1 = list()
        line_chunk2 = list()
        for line1, line2 in zip(self.openfile, self.openfile2):
            if self.gzip:
                line1 = line1.decode("utf-8")
                line2 = line2.decode("utf-8")
            line_chunk1.append(line1)
            line_chunk2.append(line2)
            if len(line_chunk1) == 4 and len(line_chunk2) == 4:
                read1 = FastqRead(line_chunk1)
                read2 = FastqRead(line_chunk2)

                # Error handling
                if not read1.header.split()[0] == read2.header.split()[0]:
                    sys.exit(
                        'INPUT ERROR: Paired reads headers does not match.\nINPUT ERROR: Read pair number:\t' +
                        '\nINPUT ERROR: ' + str(read1.header) + '\nINPUT ERROR: ' + str(read2.header) +
                        '\nINPUT ERROR: Exiting')
                line_chunk1 = list()
                line_chunk2 = list()
                yield read1, read2

    def close(self):
        """
        Closes files properly so they can be re-read if need be.
        :return: None
        """
        self.openfile.close()
        if self.filehandle2:
            self.openfile2.close()


class FastqRead(object):
    """
    Stores read as instance.
    """

    def __init__(self, fastq_as_line):
        """
        Setup function, creates read objects from lines (read = FastqRead(four_lines)), will have variables .header,
        .seq, .comment and .qual.
        :param fastq_as_line: string. Four lines (separated by newline) in fastq format.
        :return: instance. Fastq read instance.
        """
        self.header = fastq_as_line[0].strip()
        self.seq = fastq_as_line[1].strip()
        self.comment = fastq_as_line[2].strip()
        self.qual = fastq_as_line[3].strip()

    def fastq_string(self):
        """
        Makes a ready-printable string from a fastq read instance.
        :return: string.
        """
        return self.header + '\n' + self.seq + '\n' + self.comment + '\n' + self.qual + '\n'


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
    max_name_width = max(map(len, summary.keys()))
    width = value_width + max_name_width + 1

    # Header
    print("="*width, file=print_to)
    print(f"STATS SUMMARY - {name}", file=print_to)
    print("-"*width, file=print_to)

    # Print stats in columns
    for name, value in summary.items():
        if type(value) is int:
            print(f"{name:<{max_name_width}} {value:>{value_width},}", file=print_to)
        elif type(value) is float:
            print(f"{name:<{max_name_width}} {value:>{value_width}.3f}", file=print_to)
        else:
            print(f"{name:<{max_name_width}} {value:>{value_width}}", file=print_to)
    print("="*width, file=print_to)
