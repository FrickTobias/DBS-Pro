#!/bin/bash
set -xoeu pipefail

rm -rf outdir

pytest -v tests/

dbspro init testdata/reads-10k-DBS.fastq.gz outdir

pushd outdir

dbspro run

m=$(cut -f 1 data.tsv | sort | md5sum | cut -f 1 -d" ")
test $m == 60df6163dafeb484e6f8746f168b9a74