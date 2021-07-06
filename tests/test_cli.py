from pathlib import Path
import pytest
import hashlib

from xopen import xopen
from dbspro.__main__ import main as dbspro_main
from dbspro.cli.init import init
from dbspro.cli.run import run
from dbspro.cli.config import change_config

TESTDATA_READS = Path("testdata/reads-10k-DBS.fastq.gz")
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
    init(workdir, [TESTDATA_READS], ABC_SEQUENCES)


def test_init_from_csv(tmpdir):
    workdir = tmpdir / "analysis"
    sample_name = "mysample"
    sample_csv = tmpdir / "samples.csv"
    with sample_csv.open(mode="w") as f:
        print(f"{TESTDATA_READS},{sample_name}", file=f)

    init(workdir, [], ABC_SEQUENCES, sample_csv=sample_csv)
    assert (workdir / sample_name + ".fastq.gz").exists()


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

    m = md5sum(workdir / "data.tsv.gz")
    assert m == "9bef4dfad7df7c09575b1ab6423f1748"


def test_version_exit_code_zero():
    with pytest.raises(SystemExit) as e:
        dbspro_main(["--version"])
    assert e.value.code == 0
