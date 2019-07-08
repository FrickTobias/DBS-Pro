#!bin/bash
set -xoeu pipefail

rm -rf outdir
bash DBS-Pro_automation.sh testdata/reads-10k-DBS.fastq.gz outdir
m=$(cat outdir/umi-counts.txt | md5sum | cut -f 1 -d" ")
test $m == 7cc629068f91a8f635a9c1e557d404e5