#!/bin/bash
set -xoeu pipefail

rm -rf outdir

dbspro run -d outdir -f testdata/reads-10k-DBS.fastq.gz \
    ../construct-info/handles.tsv ../construct-info/ABC-sequences.fasta

m=$(cat outdir/umi-counts.txt | md5sum | cut -f 1 -d" ")
test $m == 62ce11e3ac4fb7428a6bf4bfed84f529