"""
Prints hello in a friendly manner.
"""

def main(args) -> None:
    print(args.string)

def add_arguments(parser):
    parser.add_argument("string", type=str, help="string")