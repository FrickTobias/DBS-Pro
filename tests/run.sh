#!bin/bash
set -xeu pipefail

rm -rf outdir
bash DBSpro_automation.sh testdata/reads.10k-DBS.fastq.gz outdir