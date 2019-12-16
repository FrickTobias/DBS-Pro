#!/bin/bash
set -xoeu pipefail

rm -rf outdir

pytest tests

dbspro init testdata/reads-10k-DBS.fastq.gz outdir

pushd outdir

dbspro run

m=$(cat umi_counts.tsv | sort | md5sum | cut -f 1 -d" ")
test $m == 3055459d0b4a38c070593438f50d442b