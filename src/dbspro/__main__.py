"""
DBS-Pro is a pipeline for processing DBS-Pro data.
"""
import sys
import logging
import pkgutil
import importlib
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import dbspro.cli as cli_package
from dbspro import __version__

logger = logging.getLogger(__name__)


def main(commandline_args=None) -> int:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(module)s - %(levelname)s: %(message)s",
                        datefmt='%Y-%m-%d %H:%M:%S')
    parser = ArgumentParser(description=__doc__, prog="dbspro")
    parser.add_argument("-v", "--version", action="version", version=__version__)
    subparsers = parser.add_subparsers()

    # Import each module that implements a subcommand and add a subparser for it.
    # Each subcommand is implemented as a module in the cli subpackage.
    # It needs to implement an add_arguments() and a main() function.
    modules = pkgutil.iter_modules(cli_package.__path__)
    for _, module_name, _ in modules:
        module = importlib.import_module("." + module_name, cli_package.__name__)
        help_message = module.__doc__.strip().split("\n", maxsplit=1)[0]
        subparser = subparsers.add_parser(
            module_name, help=help_message, description=module.__doc__,
            formatter_class=RawDescriptionHelpFormatter
        )
        subparser.set_defaults(module=module)
        module.add_arguments(subparser)

    args = parser.parse_args(commandline_args)
    if not hasattr(args, "module"):
        parser.error("Please provide the name of a subcommand to run")
    else:
        module = args.module
        del args.module

        # Print settings for module
        sys.stderr.write(f"SETTINGS FOR: {module.__name__.split('.')[-1]} (version: {__version__})\n")
        for object_variable, value in vars(args).items():
            sys.stderr.write(f" {object_variable}: {value}\n")

        module.main(args)

    return 0


if __name__ == "__main__":
    sys.exit(main())
