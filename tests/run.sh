#!bin/bash
set -xoeu pipefail

rm -rf outdir

dbspro run -d outdir -f testdata/reads-10k-DBS.fastq.gz \
    construct-info/handles.tsv construct-info/ABC-sequences.tsv

m=$(cat outdir/umi-counts.txt | md5sum | cut -f 1 -d" ")
test $m == d7e4d6d91069a97e2fd2cc3945222167