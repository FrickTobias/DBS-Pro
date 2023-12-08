from pathlib import Path
import pytest
import os

import pandas as pd
from dbspro.__main__ import main as dbspro_main
from dbspro.cli.init import init
from dbspro.cli.run import run
from dbspro.cli.config import run_config, load_yaml

TESTDATA_DIR = Path("testdata")
DBS_PRO_V1_DIR = TESTDATA_DIR / "dbspro_v1"
DBS_PRO_V1_SAMPLE1_READS = DBS_PRO_V1_DIR / "sample1.fastq.gz"
DBS_PRO_V1_SAMPLE2_READS = DBS_PRO_V1_DIR / "sample2.fastq.gz"
DBS_PRO_V1_SAMPLE3_READS = DBS_PRO_V1_DIR / "sample3.fastq.gz"
DBS_PRO_V1_ABC_SEQUENCES = DBS_PRO_V1_DIR / "ABC-sequences.fasta"

DBS_PRO_V3_DIR = TESTDATA_DIR / "dbspro_v3"
DBS_PRO_V3_SAMPLE_READS = DBS_PRO_V3_DIR / "sample.fastq.gz"
DBS_PRO_V3_ABC_SEQUENCES = DBS_PRO_V3_DIR / "ABCs.fasta"

EXPECTED_OUTPUT_FILES = ["report.ipynb", "data.tsv.gz", "multiqc_report.html"]


@pytest.fixture(scope="session", autouse=True)
def execute_before_any_test():
    # Set PYTHONHASHSEED to 1 to make sure that the hash of a dictionary is the same
    # Usually UMI-tools will differ in output unless this is set
    os.environ["PYTHONHASHSEED"] = "1"


def test_init(tmpdir):
    workdir = tmpdir / "analysis"
    init(workdir, [DBS_PRO_V1_SAMPLE1_READS], DBS_PRO_V1_ABC_SEQUENCES)


def test_init_from_csv(tmpdir):
    workdir = tmpdir / "analysis"
    sample_name = "mysample"
    sample_csv = tmpdir / "samples.csv"
    with sample_csv.open(mode="w") as f:
        print(f"{DBS_PRO_V1_SAMPLE1_READS},{sample_name}", file=f)

    init(workdir, [], DBS_PRO_V1_ABC_SEQUENCES, sample_csv=sample_csv)
    assert (workdir / sample_name + ".fastq.gz").exists()


def test_change_config(tmpdir):
    workdir = tmpdir / "analysis"
    init(workdir, [DBS_PRO_V1_SAMPLE1_READS], DBS_PRO_V1_ABC_SEQUENCES)
    run_config(workdir / "dbspro.yaml", [("dbs_cluster_dist", "1")])


def test_change_config_construct(tmpdir):
    workdir = tmpdir / "analysis"
    init(workdir, [DBS_PRO_V1_SAMPLE1_READS], DBS_PRO_V1_ABC_SEQUENCES)
    # Change to construct 1
    run_config(workdir / "dbspro.yaml", construct="dbspro_v1")
    configs1, _ = load_yaml(workdir / "dbspro.yaml")

    # Change to construct 3
    run_config(workdir / "dbspro.yaml", construct="dbspro_v3")
    config3, _ = load_yaml(workdir / "dbspro.yaml")
    assert configs1["h1"] != config3["h1"]


@pytest.fixture(scope="module")
def workdir_dbspro_v1(tmp_path_factory):
    """This runs the pipeline using default parameters"""
    path = tmp_path_factory.mktemp(basename="analysis-")
    workdir = path / "analysis"
    sample_csv = path / "samples.csv"
    with sample_csv.open(mode="w") as f:
        print(f"{DBS_PRO_V1_SAMPLE1_READS},Sample1", file=f)
        print(f"{DBS_PRO_V1_SAMPLE2_READS},Sample2", file=f)
        print(f"{DBS_PRO_V1_SAMPLE3_READS},Sample3", file=f)

    init(workdir, [], DBS_PRO_V1_ABC_SEQUENCES, sample_csv=sample_csv)
    run_config(workdir / "dbspro.yaml", construct="dbspro_v1")
    run(cores=4, workdir=workdir)
    return workdir


@pytest.fixture(scope="module")
def workdir_dbspro_v3(tmp_path_factory):
    """This runs the pipeline using default parameters"""
    path = tmp_path_factory.mktemp(basename="analysis-")
    workdir = path / "analysis"
    sample_csv = path / "samples.csv"
    with sample_csv.open(mode="w") as f:
        print(f"{DBS_PRO_V3_SAMPLE_READS},Sample", file=f)

    init(workdir, [], DBS_PRO_V3_ABC_SEQUENCES, sample_csv=sample_csv)
    run_config(workdir / "dbspro.yaml", construct="dbspro_v3")
    run(cores=4, workdir=workdir)
    return workdir


@pytest.mark.parametrize("file", EXPECTED_OUTPUT_FILES)
def test_output_files_exists_dbspro_v1(workdir_dbspro_v1, file):
    assert (workdir_dbspro_v1 / file).is_file()


@pytest.mark.parametrize("file", EXPECTED_OUTPUT_FILES)
def test_output_files_exists_dbspro_v3(workdir_dbspro_v3, file):
    assert (workdir_dbspro_v3 / file).is_file()


def test_output_tsv_correct_dbspro_v1(workdir_dbspro_v1):
    data = pd.read_csv(workdir_dbspro_v1 / "data.tsv.gz", sep="\t")

    assert len(set(data["Sample"])) == 3, "Should be 3 samples"
    assert len(set(data["Target"])) == 3, "Should be 3 targets"
    assert len(data) == 2499, "Should be 2495 rows"
    assert sum(data["ReadCount"]) == 52998, "Should be 52998 reads"
    assert len(set(data["Barcode"])) == 2378, "Should be 2378 unique barcodes"


def test_output_tsv_correct_dbspro_v3(workdir_dbspro_v3):
    data = pd.read_csv(workdir_dbspro_v3 / "data.tsv.gz", sep="\t")

    assert len(set(data["Sample"])) == 1, "Should be 3 samples"
    assert len(set(data["Target"])) == 7, "Should be 3 targets"
    assert len(data) == 3374, "Should be 3374 rows"
    assert sum(data["ReadCount"]) == 16482, "Should be 16482 reads"
    assert len(set(data["Barcode"])) == 2287, "Should be 2378 unique barcodes"


def test_version_exit_code_zero():
    with pytest.raises(SystemExit) as e:
        dbspro_main(["--version"])
    assert e.value.code == 0
