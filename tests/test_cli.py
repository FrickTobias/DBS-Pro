from pathlib import Path
import pytest
import hashlib
import os

from xopen import xopen
from dbspro.__main__ import main as dbspro_main
from dbspro.cli.init import init
from dbspro.cli.run import run
from dbspro.cli.config import change_config

TESTDATA_SAMPLE1_READS = Path("testdata/sample1.fastq.gz")
TESTDATA_SAMPLE2_READS = Path("testdata/sample2.fastq.gz")
TESTDATA_SAMPLE3_READS = Path("testdata/sample3.fastq.gz")
ABC_SEQUENCES = Path("testdata/ABC-sequences.fasta")


def md5sum(filename):
    md5_hash = hashlib.md5()
    with xopen(filename, mode="rb") as f:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            md5_hash.update(byte_block)
    return md5_hash.hexdigest()


def test_init(tmpdir):
    workdir = tmpdir / "analysis"
    init(workdir, [TESTDATA_SAMPLE1_READS], ABC_SEQUENCES)


def test_init_from_csv(tmpdir):
    workdir = tmpdir / "analysis"
    sample_name = "mysample"
    sample_csv = tmpdir / "samples.csv"
    with sample_csv.open(mode="w") as f:
        print(f"{TESTDATA_SAMPLE1_READS},{sample_name}", file=f)

    init(workdir, [], ABC_SEQUENCES, sample_csv=sample_csv)
    assert (workdir / sample_name + ".fastq.gz").exists()


def test_change_config(tmpdir):
    workdir = tmpdir / "analysis"
    init(workdir, [TESTDATA_SAMPLE1_READS], ABC_SEQUENCES)
    change_config(workdir / "dbspro.yaml", [("dbs_cluster_dist", "1")])


@pytest.fixture(scope="module")
def workdir(tmp_path_factory):
    """This runs the pipeline using default parameters"""
    path = tmp_path_factory.mktemp(basename="analysis-")
    workdir = path / "analysis"
    sample_csv = path / "samples.csv"
    with sample_csv.open(mode="w") as f:
        print(f"{TESTDATA_SAMPLE1_READS},Sample1", file=f)
        print(f"{TESTDATA_SAMPLE2_READS},Sample2", file=f)
        print(f"{TESTDATA_SAMPLE3_READS},Sample3", file=f)

    init(workdir, [], ABC_SEQUENCES, sample_csv=sample_csv)
    run(targets=[], workdir=workdir)
    return workdir


@pytest.mark.parametrize("file", ["report.ipynb", "data.tsv.gz"])
def test_output_files_exists(workdir, file):
    assert (workdir / file).is_file()


def test_output_tsv_md5sum(workdir):
    if "PYTHONHASHSEED" not in os.environ or os.environ["PYTHONHASHSEED"] != "1":
        pytest.skip("Environment variable 'PYTHONHASHSEED' not '1'.")
    else:
        m = md5sum(workdir / "data.tsv.gz")
        assert m == "b908621c64f067a3e0c1a986bee0fdcb"


def test_version_exit_code_zero():
    with pytest.raises(SystemExit) as e:
        dbspro_main(["--version"])
    assert e.value.code == 0
