"""
Set DBS-Pro to work on custom read structures
"""

import logging
import os
import sys

logger = logging.getLogger(__name__)


def main(args):

    if args.print_structure:
        print("h1-DBS-h2-ABC-UMI-h3\n"
              "DBS length: 20\n"
              "UMI length: 6")
        sys.exit()

    outline = str()
    unchanged = True
    infile = file_name_fetcher(construct=args.construct)
    with open(infile, "r") as openin:
        for line in openin:
            key, value = line.split()

            # Only prints current sequences
            if args.show_sequences:
                logger.info(key + ' ' + value)

            # Change the construct specified to new seq
            elif args.construct == key:
                logger.info("Changing " + args.construct)
                logger.info("Prev sequence: " + value)
                logger.info("New sequence: " + args.sequence)
                if value == args.sequence:
                   logger.warning("Previous and new sequence is identical")

                value = args.sequence
                unchanged = False

            outline += (key + "\t" + value + "\n")

    if unchanged:
        logger.error(f"Invalid construct name, file unchanged ({args.construct})")
        sys.exit()
    else:
        with open(infile, 'w') as openout:
            openout.write(outline)


def file_name_fetcher(construct):
    """
    Fetches the path to the file which should be altered depending on the first letter of the
    :return: file_path
    """

    file_dict = {
        "h":"/construct-info/handles.tsv",
        "A":"/construct-info/ABC-sequences.tsv"
    }

    dbs_pro_folder = os.path.realpath(__file__).rsplit("/",4)[0]
    try: filename = dbs_pro_folder + file_dict[construct[0]]
    except KeyError:
        logger.error(f"Invalid construct name, file unchanged ({construct})")
        sys.exit()
    return filename

def add_arguments(parser):
    parser.add_argument("construct", help="Part of the construct to be changed. Possible values: h1, h2, h3, ABC1, "
                                          "ABC2, ABC3.")
    parser.add_argument("sequence", help="Sequence to be changed to, should be written 5' to 3'")
    parser.add_argument("-p", "--print_structure", action="store_true", help="Prints the assumed read structure and "
                                                                             "exits")
    parser.add_argument("-s","--show_sequences", action="store_true", help="Show all sequences currently set for the "
                                                                           "handles or ABC, depending on the <construct> "
                                                                           "argument")