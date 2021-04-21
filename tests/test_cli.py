from pathlib import Path

from dbspro.cli.init import init
from dbspro.cli.run import run
from dbspro.cli.config import change_config

TESTDATA_READS = Path("testdata/reads-10k-DBS.fastq.gz")
ABC_SEQUENCES = Path("testdata/ABC-sequences.fasta")
DBSPRO_LOG_FILE = Path("tests/data/example_dbspro_logfile.log")
CUTADAPT_LOG_FILE = Path("tests/data/example_cutadapt_logfile.log")


def test_init(tmpdir):
    workdir = tmpdir / "analysis"
    init(workdir, [TESTDATA_READS], ABC_SEQUENCES)


def test_change_config(tmpdir):
    workdir = tmpdir / "analysis"
    init(workdir, [TESTDATA_READS], ABC_SEQUENCES)
    change_config(workdir / "dbspro.yaml", [("dbs_cluster_dist", "1")])


def test_output_files(tmpdir, targets=["report.ipynb", "data.tsv.gz"]):
    workdir = tmpdir / "analysis"
    init(workdir, [TESTDATA_READS], ABC_SEQUENCES)
    run(targets=targets, workdir=workdir)
    for target in targets:
        assert Path(f"{workdir}/{target}").is_file()
