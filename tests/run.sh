#!/bin/bash
set -xoeu pipefail

rm -rf outdir

dbspro init testdata/reads-10k-DBS.fastq.gz outdir

pushd outdir

dbspro run

m=$(cat outdir/umi-counts.txt | sort | md5sum | cut -f 1 -d" ")
test $m == 3c31f9634af4238cb0037170104e6d5e