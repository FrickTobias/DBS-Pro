#! /usr/bin/env python2

def main():

    #
    # Imports & globals
    #
    global args, summaryInstance, sys, time, script_name
    import sys, time

    script_name = sys.argv[0].split('/')[-1].split('.')[0]

    #
    # Argument parsing
    #
    argumentsInstance = readArgs()

    #
    # Process data
    #

    # Barcode processing
    report_progress("Starting analysis")
    report_progress("Saving DBS information to RAM")
    progress = ProgressReporter("Lines read", 1000000)
    generator = FileReader(args.dbs)
    bc_dict = dict()
    for read in generator.fastqReader():
        bc_dict[read.header] = read.seq
        progress.update()
    report_progress("Finished processing DBS sequences")

    # Counting UMI:s found in the different ABC:s for all barcodes.
    report_progress("Calculating stats")
    result_dict = dict()
    abc_list = [args.umi_1, args.umi_2, args.umi_3]
    umi_without_proper_bc = int()
    for current_abc in abc_list:
        report_progress("Reading file\t" + str(current_abc))
        generator = FileReader(current_abc)
        progress = ProgressReporter("Lines read", 1000000)

        # Loop over reads in file, where read.seq = umi
        for read in generator.fastqReader():

            # Try find UMI
            try: bc = bc_dict[read.header]
            except KeyError:
                umi_without_proper_bc += 1
                progress.update()
                continue

            # If not dbs in result dict, add it and give it a dictionary for every abc
            if not bc in result_dict:
                result_dict[bc] = dict()
                for abc in abc_list:
                    result_dict[bc][abc] = dict()

            # Add +1 to corresponding UMI sequence
            if not read.seq in result_dict[bc][abc]:
                result_dict[bc][current_abc][read.seq] = int()
            result_dict[bc][current_abc][read.seq] += 1

            progress.update()
        report_progress("Finished reading file\t" + str(current_abc))

    # Barcode-globbing umi/read counter for all ABC:s
    abc_counter_umi = dict()
    abc_counter_read = dict()
    for abc in abc_list:
        abc_counter_umi[abc] = [0]
        abc_counter_read[abc] = [0]

    # Output file writing and
    with open(args.output, 'w') as openout:
        for bc in result_dict.keys():
            out_string = str()
            for abc in abc_list:
                # Prepping outstring: umi_count(abc1) + \t + umi_count(abc2) + \t + umi_count(abc3) \n
                out_string += str(len(result_dict[bc][abc])) + '\t'
                # Add number of UMI:s
                if sum(result_dict[bc][abc].values()) >= args.filter:
                    abc_counter_umi[abc].append(len(result_dict[bc][abc].keys()))
                    # Add number of reads
                    abc_counter_read[abc].append(sum(result_dict[bc][abc].values()))
                    # If not enough reads, remove entry
                #else:
                #    del result_dict[bc][abc]

            openout.write(out_string + '\n')

    # Reporting stats to terminal
    print()
    report_progress("Tot DBS count:\t" + str(len(result_dict.keys())))
    for abc in abc_list:
        print()
        report_progress("\t" + abc + " tot umi count:\t" + "{:,}".format(sum(abc_counter_umi[abc])))
        report_progress("\t" + abc + " N50 umi per dbs:\t" + "{:,}".format(n50_counter(abc_counter_umi[abc])))
        report_progress("\t" + abc + " tot read count:\t" + "{:,}".format(sum(abc_counter_read[abc])))
        report_progress("\t" + abc + " N50 read per dbs:\t" + "{:,}".format(n50_counter(abc_counter_read[abc])))
    print()

    # Plotting
    report_progress("Prepping data for plot")
    read_dict_for_plotting, umi_dict_for_plotting = format_data_for_plotting(result_dict)
    plot_density_correlation_matrix(args.read_plot,read_dict_for_plotting)
    plot_density_correlation_matrix(args.umi_plot,umi_dict_for_plotting)
    report_progress("Finished")

def dict_clearer(dictionary):
    """
    Takes a dictionary and removed any keys which does not have any values
    :param dictionary:
    :return:
    """

    print(dictionary)

    for key, val in dictionary.copy().items():
        if len(val) == 0:
            del dictionary[key]

    return dictionary

def n50_counter(input_list):
    """
    Calculates N50 for a given list
    :param list: list with numbers (list)
    :return: N50 (same type as elements of list)
    """
    input_list.sort()
    half_tot = sum(input_list)/2

    current_count = 0
    for num in input_list:
        current_count += num
        if current_count >= half_tot:
            return num

def format_data_for_plotting(result_dict):
    """
    Divides result dict to plot-ready dicts
    :param result_dict: dict[dbs][abc][umi] = read_count
    :return:    dict[dbs][abc] = read_count
                dict[dbs][abc] = umi_count
    """
    read_dict_for_plotting = dict()
    umi_dict_for_plotting = dict()
    for dbs in result_dict.keys():
        read_dict_for_plotting[dbs] = dict()
        umi_dict_for_plotting[dbs] = dict()
        for abc in result_dict[dbs].keys():
            read_dict_for_plotting[dbs][abc] = sum(result_dict[dbs][abc].values())
            umi_dict_for_plotting[dbs][abc] = len(result_dict[dbs][abc].keys())

    return read_dict_for_plotting, umi_dict_for_plotting

def plot_density_correlation_matrix(name, result_dict):
    """

    :param x:
    :param y:
    :param z:
    :return:
    """


    # With two column indices, values same
    # as dictionary keys
    import pandas as pd

    df_list = list()
    for dbs_dict in result_dict.values():
        df_list.append(dbs_dict)

    df = pd.DataFrame(df_list, columns=[args.umi_1, args.umi_2, args.umi_3])

    # Make plot imports
    import seaborn as sns
    import matplotlib.pyplot as plt

    # Basic correlogram
    #g = sns.pairplot(df, diag_kws=dict(bins=50))
    g = sns.pairplot(df, diag_kind="kde", diag_kws=dict(shade=True, bw=.05, vertical=False))
    #g.map_diag(sns.kdeplot, lw = 3, legend = False)#plt.hist, bins=200)
    #g.map_diag(plt.hist)

    for x in range(3):
        for y in range(3):
            g.axes[x,y].set_xlim((0, 50))
            g.axes[x,y].set_ylim((0, 50))
    report_progress("Plotting " + name)
    import os
    my_path = os.path.abspath(__file__)
    plt.savefig(name)

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

class readArgs(object):
    """
    Reads arguments and handles basic error handling like python version control etc.
    """

    def __init__(self):

        readArgs.parse(self)
        readArgs.pythonVersion(self)

    def parse(self):

        #
        # Imports & globals
        #
        import argparse
        global args

        parser = argparse.ArgumentParser(description=__doc__)

        # Arguments
        parser.add_argument("dbs", help="Reads with only DBS seq in fastq format.")
        parser.add_argument("umi_1", help="Reads with only UMI seq (unique molecular identifier) file for ABC (antibody barcode) 1 in fastq format")
        parser.add_argument("umi_2", help="Reads with only UMI seq (unique molecular identifier) file for ABC (antibody barcode) 2 in fastq format")
        parser.add_argument("umi_3", help="Reads with only UMI seq (unique molecular identifier) file for ABC (antibody barcode) 3 in fastq format")
        parser.add_argument("output", help="output file")
        parser.add_argument("read_plot", help="Filename for output reads/DBS pair plot (will be .png)")
        parser.add_argument("umi_plot", help="Filename for output UMI:s/DBS pair plot (will be .png)")
        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3.")
        parser.add_argument("-f", "--filter", type=int, default=0, help="Number of minimum reads required for an ABC "
                                                                        "to be included in output. DEFAULT: 0")

        args = parser.parse_args()

    def pythonVersion(self):
        """ Makes sure the user is running python 3."""

        #
        # Version control
        #
        import sys
        if sys.version_info.major == 3:
            pass
        else:
            sys.stderr.write('\nWARNING: you are running python ' + str(
                sys.version_info.major) + ', this script is written for python 3.')
            if not args.force_run:
                sys.stderr.write('\nAborting analysis. Use -F (--Force) to run anyway.\n')
                sys.exit()
            else:
                sys.stderr.write('\nForcing run. This might yield inaccurate results.\n')

#class Summary(object):
#
#    def __init__(self):
#
#        self.variable = int()
#
#    def writeToStdErr(self):
#        """
#        Writes all object variables to stdout.
#        """
#
#        for objectVariable, value in vars(self).items():
#            sys.stderr.write('\n\n' + str(objectVariable) + '\n' + str(value))
#        sys.stderr.write('\n')
#
#    def writeLog(self):
#        """
#        Writes all object variables to a log file (outfile.log)
#        """
#
#        self.log = args.outfile + '.log'
#        import time
#        with open(self.log, 'w') as openout:
#            openout.write(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
#            for objectVariable, value in vars(self).items():
#                openout.write('\n\n'+str(objectVariable) + '\n' + str(value))

if __name__=="__main__": main()