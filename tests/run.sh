#!bin/bash
set -xoeu pipefail

rm -rf outdir

dbspro run -d outdir -f testdata/reads-10k-DBS.fastq.gz

m=$(cat outdir/umi-counts.txt | md5sum | cut -f 1 -d" ")
test $m == 3c086da27caec5dd08d0da695b81b067