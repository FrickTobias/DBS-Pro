#!bin/bash
set -xeu pipefail

rm -rf outdir
bash DBSpro_automation.sh DBSpro-testdata-0.1/2-PrEST02_S11_L001_R1_001.fastq.gz outdir