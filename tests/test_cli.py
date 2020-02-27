from pathlib import Path
from io import StringIO
import pytest

from dbspro.cli.init import init
from dbspro.cli.run import run
from dbspro.cli.config import change_config
from dbspro.cli.summary import Summary, parse_dbspro, parse_cutadapt

TESTDATA_READS = Path("testdata/reads-10k-DBS.fastq.gz")
ABC_SEQUENCES = Path("testdata/ABC-sequences.fasta")
DBSPRO_LOG_FILE = Path("tests/data/example_dbspro_logfile.log")
CUTADAPT_LOG_FILE = Path("tests/data/example_cutadapt_logfile.log")


def test_init(tmpdir):
    workdir = tmpdir / "analysis"
    init(workdir, TESTDATA_READS, ABC_SEQUENCES)


def test_change_config(tmpdir):
    workdir = tmpdir / "analysis"
    init(workdir, TESTDATA_READS, ABC_SEQUENCES)
    change_config(workdir / "dbspro.yaml", [("dbs_cluster_dist", "1")])


def test_output_files(tmpdir, targets=["report.ipynb", "data.tsv"]):
    workdir = tmpdir / "analysis"
    init(workdir, TESTDATA_READS, ABC_SEQUENCES)
    run(targets=targets, workdir=workdir)
    for target in targets:
        assert Path(f"{workdir}/{target}").is_file()


def test_summary_parse_dbspro():
    summary = Summary()
    parse_dbspro(DBSPRO_LOG_FILE, summary)

    true_values = [
        ("Reads kept", "9592"),
        ("DBS clusters linked to ABC", "2932"),
        ("Total UMIs", "3312"),
        ("Total clustered UMIs", "5478")
    ]

    for filename, values in summary.data.items():
        assert filename == DBSPRO_LOG_FILE.stem
        for (parameter, value), (true_param, true_value) in zip(values.items(), true_values):
            assert parameter == true_param
            assert value == true_value


def test_summary_parse_cutadapt():
    summary = Summary()
    parse_cutadapt(CUTADAPT_LOG_FILE, summary)

    true_values = [
        ("Total reads processed", "100290"),
        ("Reads with adapters", "100290"),
        ("Reads with adapters (%)", "100.0"),
        ("Reads that were too short", "95"),
        ("Reads that were too short (%)", "0.1"),
        ("Reads that were too long", "13"),
        ("Reads that were too long (%)", "0.0"),
        ("Reads written (passing filters)", "100182"),
        ("Reads written (passing filters) (%)", "99.9")
    ]

    for filename, values in summary.data.items():
        assert filename == CUTADAPT_LOG_FILE.stem
        for (parameter, value), (true_param, true_value) in zip(values.items(), true_values):
            assert parameter == true_param
            assert value == true_value
