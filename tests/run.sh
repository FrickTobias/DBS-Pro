#!bin/bash
set -xoeu pipefail

rm -rf outdir
bash DBSpro_automation.sh testdata/reads.10k-DBS.fastq.gz outdir
m=$(cat outdir/umi-counts.txt | md5sum | cut -f 1 -d" ")
test $m == 3c086da27caec5dd08d0da695b81b067