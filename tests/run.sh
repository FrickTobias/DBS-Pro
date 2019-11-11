#!/bin/bash
set -xoeu pipefail

rm -rf outdir

dbspro init testdata/reads-10k-DBS.fastq.gz outdir

pushd outdir

dbspro run

m=$(cat outdir/umi-counts.txt | md5sum | cut -f 1 -d" ")
test $m == 62ce11e3ac4fb7428a6bf4bfed84f529