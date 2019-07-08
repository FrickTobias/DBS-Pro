"""
Set DBS-Pro to work on custom read structures
"""

import logging
import os
import sys

logger = logging.getLogger(__name__)


def main(args):

    # Print constrict info
    if args.print_structure:
        print("h1-DBS-h2-ABC-UMI-h3\n"
              "DBS length: 20\n"
              "UMI length: 6")
        sys.exit()

    # Read existing file
    unchanged = True
    infile = file_name_fetcher(construct=args.construct)
    outdict = dict()
    with open(infile, "r") as openin:
        for line in openin:
            key, value = line.split()
            if args.show_sequences:
                logger.info(f"{key} {value}")
            else:
                outdict[key] = value

    # Modify file
    if not args.show_sequences:

        # Modify dictionary
        if args.construct in outdict:
            # Remove
            if args.remove:
                logger.info(f"Removing {args.construct} from {infile}")
                del outdict[args.construct]
                unchanged = False
            # Change
            elif outdict[args.construct] == args.sequence:
                pass
            else:
                logger.info(f"Changing {args.construct} to {args.sequence}")
                outdict[args.construct] = args.sequence
                unchanged = False
        # Add
        elif args.add:
            logger.info(f"Adding {args.construct} {args.sequence} to {infile}")
            outdict[args.construct] = args.sequence
            unchanged = False

        # Write output
        if unchanged:logger.error(f"File has not been changed, please review your arguments. (construct: "
                                  f"{args.construct}, seuqence: {args.sequence})")
        else:
            with open(infile, 'w') as openout:
                openout.write("Antibody-target\tBarcode-sequence\n")
                for key, value in outdict.items():
                    openout.write(f"{key}\t{value}\n")


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
    parser.add_argument("-a","--add", action="store_true", help="Add ABC sequence instead of change it. Will fail if "
                                                                "used with an already existing construct name.")
    parser.add_argument("-r","--remove", action="store_true", help="Remove ABC sequence instead of change it. Will fail "
                                                                   "if construct name is not found.")