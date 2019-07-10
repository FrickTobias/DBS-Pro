#!bin/bash
set -xoeu pipefail

rm -rf outdir

dbspro run -d outdir -f testdata/reads-10k-DBS.fastq.gz

m=$(cat outdir/umi-counts.txt | md5sum | cut -f 1 -d" ")
test $m == ae75d387667946c51f734d4de9f843b5