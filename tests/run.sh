#!/bin/bash
set -xoeu pipefail

rm -rf outdir

pytest -v tests/

dbspro init testdata/reads-10k-DBS.fastq.gz outdir --abc testdata/ABC-sequences.fasta

pushd outdir

dbspro run

m=$(cut -f 1 data.tsv | sort | md5sum | cut -f 1 -d" ")
test $m == f2703c12e53b56d32d3cdaee017728e1