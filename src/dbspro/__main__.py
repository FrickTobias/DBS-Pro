"""
A pipeline for processing DBS-Pro data.

Workflow:

1. Initiate a new analysis:

    dbspro init -h

2. Change the config file to your liking:

    dbspro config -h

3. Run the pipeline:

    dbspro run -h

For more information about the pipeline, see the documentation at
https://github.com/FrickTobias/DBS-Pro

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
    parser = ArgumentParser(description=__doc__, prog="dbspro", formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--version", action="version", version=__version__)
    parser.add_argument("--profile", action="store_true", default=False,
                        help="Save profiling info to dbspro_<subcommand>.prof")
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

    args, extra_args = parser.parse_known_args(commandline_args)
    if not hasattr(args, "module"):
        parser.error("Please provide the name of a subcommand to run")

    module = args.module
    subcommand = module.main
    del args.module
    profile = args.profile
    del args.profile

    module_name = module.__name__.split('.')[-1]

    if extra_args:
        # Re-parse extra arguments if module is not "run" to raise the expected error
        if module_name != "run":
            parser.parse_args(extra_args)

        # Add extra arguments to args.snakemake_args for module "run"
        args.snakemake_args = extra_args

    # Print settings for module
    sys.stderr.write(f"SETTINGS FOR: {module_name} (version: {__version__})\n")
    for object_variable, value in vars(args).items():
        sys.stderr.write(f" {object_variable}: {value}\n")

    if profile:
        import cProfile
        profile_file = f'dbspro_{module_name}.prof'
        cProfile.runctx("subcommand(args)", globals(), dict(subcommand=subcommand, args=args),
                        filename=profile_file)
        logger.info(f"Writing profiling stats to '{profile_file}'.")
    else:
        subcommand(args)

    return 0


if __name__ == "__main__":
    sys.exit(main())
