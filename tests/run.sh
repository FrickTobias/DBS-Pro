#!/bin/bash
set -xoeu pipefail

rm -rf outdir

pytest -v tests/

dbspro init testdata/reads-10k-DBS.fastq.gz outdir --abc testdata/ABC-sequences.fasta

pushd outdir

dbspro run

m=$(cut -f 1 data.tsv | sort | md5sum | cut -f 1 -d" ")
test $m == 590a6a3d182f88cad7f9404f8b7445db