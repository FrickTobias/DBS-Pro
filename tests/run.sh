#!/bin/bash
set -xoeu pipefail

if [[ $(uname) == Darwin ]]; then
  md5="md5 -r"
else
  md5=md5sum
fi

rm -rf outdir

pytest -v tests/

dbspro init outdir testdata/reads-10k-DBS.fastq.gz --abc testdata/ABC-sequences.fasta

pushd outdir

dbspro run

m=$(gunzip -c data.tsv.gz | cut -f 1 | sort | $md5 | cut -f 1 -d" ")
test $m == 590a6a3d182f88cad7f9404f8b7445db