"""
Combines starcode output files into raw fastq files creating error corrected fastq files
"""
def main(args):

    global sys, time, script_name, gzip
    import sys, time, gzip

    script_name = sys.argv[0].split('/')[-1].split('.')[0]

    report_progress("Starting analysis")
    report_progress("Processing file " + args.err_corr)
    progress = ProgressReporter("Lines parsed", 1000000)
    generator = FileReader(args.err_corr)
    err_corr = dict()
    for line in generator.fileReader():

        try: cluster_seq, num_reads, raw_seqs_list = line.split()
        except ValueError:
            print("WARNING! Non-default starcode output line:\t'" + line + "'")
        for raw_seq in raw_seqs_list.split(","):
            if not raw_seq in err_corr:
                err_corr[raw_seq] = cluster_seq
        progress.update()
    generator.close()
    report_progress("Error corrected sequenced parsed.")

    report_progress("Correcting sequences and writing to output file.")
    progress = ProgressReporter("Lines parsed", 1000000)
    generator = FileReader(args.raw_fastq)
    no_err_corr_seq = int()
    tot_reads = int()
    corr_seqs = int()
    with open(args.corr_fastq, 'w') as openout:
        for read in generator.fastqReader():

            tot_reads += 1
            if read.seq in err_corr:
                cluster_seq = err_corr[read.seq]
                read.seq = cluster_seq
                openout.write(read.fastq_string())
                corr_seqs += 1
            else:
                no_err_corr_seq += 1
            progress.update()
    generator.close()

    report_progress("Reads total:\t" + str(tot_reads))
    report_progress("Reads corrected:\t" + str(corr_seqs))
    report_progress("Reads without corrected seq:\t" + str(no_err_corr_seq))
    report_progress("Finished")

def report_progress(string):
    """
    Progress report function, writes string (<time_stamp> <tab> <string> <newline>) to std err.
    :param string: String to be written to terminal
    :return: None
    """
    sys.stderr.write(script_name.upper() + ':\t' + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + '\t' + string + '\n')

class ProgressReporter(object):
    """
    Progress reporter object, writes updates to std err for every <report_step> iteration. Used for iterations of unknown
    lengths.
    It is setup by first creating an instance of the object (progress = ProgressReporter(name_of_process='file_reading',
    report_step=1000)) and then using updating the progress every iteration (progress.update()). In the example given
    this would print <time_stamp> <tab> <'file_reading'> <tab> <iterations_made> every 1000th iteration.
    """

    def __init__(self, name_of_process, report_step):
        """
        Setup function for iterations of unknown lengths. Typically initial file readings.
        :param name_of_process: string. Name of the iterations being made
        :param report_step: integer. Will write update to terminal every <integer> interations.
        :return: None
        """

        self.name = name_of_process
        self.report_step = report_step
        self.position = int()
        self.next_limit = report_step

    def update(self):
        """
        Update function for object. Run once every iteration, does not require arguments.
        :return: None
        """
        self.position += 1
        if self.position >= self.next_limit:
            report_progress(self.name + '\t' + "{:,}".format(self.position))
            self.next_limit += self.report_step

class ProgressBar(object):
    """
    Progress reporter object, writes updates to std err in for of a progress bar. Used for iterations of known lenghts.
    It is setup by first creating an instance of the object (progressBar = ProgressBar(name='processing', min=0,
    max=100000, step=1)) and then using the updating function every iteration (progressBar.update()). This creates a
    progress bar in the terminal which updates for every 2% reached of the process.
    """
    def __init__(self, name, min, max, step):
        """
        Setup function for iterations of known lengths. Typically used for processing after having read a file.
        :param name: string. Name printed above progress bar
        :param min: integer. Starting point of process, most often 0.
        :param max: integer. Last iteration of process.
        :param step: integer. How many iterations have been made for every time the update() function is called.
        :return: None
        """

        # Variables
        self.min = min
        self.max = max
        self.current_position = min
        self.step = step

        # Metadata
        self.two_percent = (self.max-self.min)/50
        self.current_percentage = self.two_percent

        # If two percent, equivalent of one '#', is less than one step length increase the number of # written each step
        if self.two_percent < self.step and not self.max==2:
            self.progress_length = int(50/(self.max-2))
            self.progress_string = '#' * self.progress_length
        elif self.max == 2:
            self.progress_string = '#' * 25
        else:
            self.progress_string = '#'

        # Printing progress bar tot length
        report_progress(str(name))
        sys.stderr.write('\n|------------------------------------------------|\n')

    def update(self):
        """
        Update function for object. Run once every <step> iteration, does not require arguments.
        :return: None
        """
        # If progress is over 2%, write '#' to stdout
        self.current_position += self.step
        if self.current_percentage < self.current_position:
            sys.stderr.write(self.progress_string)
            sys.stderr.flush()
            time.sleep(0.001)
            self.current_percentage += self.two_percent

    def terminate(self):
        """
        Termination function. Writes newline to std err, typically used directly after iteration loop is complete.
        :return: None
        """
        sys.stderr.write('\n')

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
            import sys
            self.openfile = sys.stdin
        # Open files as zipped or not not (depending on if they end with .gz)
        elif self.filehandle[-3:] == '.gz':
            report_progress('File detected as gzipped, unzipping when reading')
            self.openfile = gzip.open(self.filehandle, 'r')
            self.gzip = True
        else:
            self.openfile = open(self.filehandle, 'r')

        # Paired end preparation
        self.filehandle2 = filehandle2
        if self.filehandle2:

            # Open files as zipped or not not (depending on if they end with .gz)
            if self.filehandle2[-3:] == '.gz':
                report_progress('File detected as gzipped, unzipping when reading')

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
                    sys.exit('INPUT ERROR: Paired reads headers does not match.\nINPUT ERROR: Read pair number:\t'+str(progress.position+1)+'\nINPUT ERROR: '+str(read1.header)+'\nINPUT ERROR: '+str(read2.header)+'\nINPUT ERROR: Exiting')
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

def add_arguments(parser):
    parser.add_argument("raw_fastq", help="Fastq file with raw sequences.")
    parser.add_argument("err_corr", help="Starcode default output with error corrected sequences.")
    parser.add_argument("corr_fastq", help="Output file in fastq with error corrected sequences.")

    # Options
    parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                       "Not recommended due to different function "
                                                                       "names in python 2 and 3.")