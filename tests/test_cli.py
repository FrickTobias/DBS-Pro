

from dbspro.cli.init import init
from dbspro.cli.run import run

TESTDATA_READS = testdata/reads-10k-DBS.fastq.gz outdir

def test_init(tmpdir):
    init(tmpdir / "analysis", TESTDATA_READS)

def test_output_files(tmpdir):
    init(tmpdir / "analysis", TESTDATA_READS)
    copy_config(("tests/test_config.yaml",
                 workdir / "blr.yaml"
                )
    run(targets)
    for target in targets
        assert Path(workdir / analysis / target).is_file()