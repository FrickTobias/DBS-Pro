#!bin/bash
set -xoeu pipefail

rm -rf outdir
bash DBS-Pro_automation.sh testdata/reads-10k-DBS.fastq.gz outdir
m=$(cat outdir/umi-counts.txt | md5sum | cut -f 1 -d" ")
test $m == ae75d387667946c51f734d4de9f843b5