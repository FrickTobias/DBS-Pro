#!/bin/bash
set -xoeu pipefail

rm -rf outdir

pytest -v tests/

dbspro init testdata/reads-10k-DBS.fastq.gz outdir

pushd outdir

dbspro run

m=$(cut -f 1,2,4 data.tsv | sort | md5sum | cut -f 1 -d" ")
test $m == e10ac1a02e57ee767619feb59ae32b21