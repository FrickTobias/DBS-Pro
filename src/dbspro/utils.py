#! /usr/bin/python3

import sys
import gzip
import logging


logger = logging.getLogger(__name__)


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
                    sys.exit('INPUT ERROR: Paired reads headers does not match.\nINPUT ERROR: Read pair number:\t'+'\nINPUT ERROR: '+str(read1.header)+'\nINPUT ERROR: '+str(read2.header)+'\nINPUT ERROR: Exiting')
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
        return self.header + '\n' + self.seq  + '\n' + self.comment  + '\n' + self.qual + '\n'
