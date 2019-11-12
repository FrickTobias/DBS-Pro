import pandas as pd
import dnaio
import os

from dbspro.utils import get_abcs

# Read sample and handles files.
abc = get_abcs("ABC-sequences.fasta")
configfile: "dbspro.yaml"

# Get required values
abc_len = list(map(len, abc['Sequence']))[0] - 1    # Assumes same length, remove one as anchored sequences
abc_umi_len = abc_len + config["UMI_len"]
dbs = "N"*config["DBS_len"]


# Define final targets for pipeline. Currently they are the output of rule 'analyze'

rule all:
    input: 'report.html'

# Cutadapt trimming

"Extract DBS and trim handle between DBS and ABC."
rule extract_dbs:
    output:
        reads="dbs-raw.fastq.gz"
    input:
        reads="reads.fastq.gz"
    log: "log_files/cutadapt-extract-dbs.log"
    threads: 20
    shell:
        "cutadapt"
        " -g ^{config[h1]}...{config[h2]}"
        " --discard-untrimmed"
        " -e 0.2"
        " -j {threads}"
        " -o {output.reads}"
        " {input.reads}"
        " > {log}"


"Extract ABC and UMI."
rule extract_abc_umi:
    output:
        reads="trimmed-abc.fastq.gz"
    input:
        reads="reads.fastq.gz"
    log: "log_files/cutadapt-extract-abc-umi.log"
    threads: 20
    shell:
        "cutadapt"
        " -g ^{config[h1]}{dbs}{config[h2]}...{config[h3]}"
        " --discard-untrimmed"
        " -e 0.2"
        " -m {abc_umi_len}"
        " -M {abc_umi_len}"
        " -j {threads}"
        " -o {output.reads}"
        " {input.reads}"
        " > {log}"

## ABCs

"Demultiplexes ABC sequnces and trims it to give ABC-specific UMI fastq files."
rule demultiplex_abc:
    output:
        reads=touch(expand("ABCs/{name}-UMI-raw.fastq.gz", name=abc['Target']))
    input:
        reads="trimmed-abc.fastq.gz"
    log: "log_files/cutadapt-id-abc.log"
    shell:
        "cutadapt"
        " -g file:ABC-sequences.fasta"
        " --no-indels"
        " -e 0.2"
        " -o ABCs/{{name}}-UMI-raw.fastq.gz"
        " {input.reads}"
        " > {log}"

# Starcode clustering

"Cluster DBS sequence using starcode"
rule dbs_cluster:
    output:
        clusters="dbs-clusters.txt"
    input:
        reads="dbs-raw.fastq.gz"
    log: "log_files/starcode-dbs-cluster.log"
    threads: 20
    shell:
        "pigz -cd {input.reads} | starcode"
        " --print-clusters"
        " -t {threads}"
        " -d {config[dbs_cluster_dist]}"
        " -o {output.clusters} 2> {log}"


"Cluster ABC sequence using starcode. If the input file is empty (no ABC sequence found)"
" a empty output file will also be created."
rule abc_cluster:
    output:
        reads="ABCs/{sample}-UMI-corrected.fasta"
    input:
        abc_reads="ABCs/{sample}-UMI-raw.fastq.gz",
        dbs_corrected="dbs-corrected.fasta"
    log: "ABCs/log_files/splitcluster-{sample}.log"
    shell:
        "dbspro splitcluster"
        " {input.dbs_corrected}"
        " {input.abc_reads}"
        " -o {output.reads}"
        " -t {config[abc_cluster_dist]}"
        " -l {config[UMI_len]} 2> {log}"


# DBS-Pro

"Combine cluster results with original files to error correct them."
rule error_correct:
    output:
        reads="{corr_file}-corrected.fasta"
    input:
        reads="{corr_file}-raw.fastq.gz",
        clusters="{corr_file}-clusters.txt"
    log: "log_files/error-correct-{corr_file}.log"
    threads: 20
    shell:
        "dbspro correctfastq"
        " {input.reads}"
        " {input.clusters}"
        " {output.reads}"

"Analyzes all result files"
rule analyze:
    output:
        umi_counts="umi_counts.tsv",
        read_counts="read_counts.tsv"
    input:
        dbs_fasta="dbs-corrected.fasta",
        abc_fastas=expand("ABCs/{abc}-UMI-corrected.fasta", abc=abc['Target'])
    log: "log_files/analyze.log"
    threads: 20
    shell:
        "dbspro analyze"
        " -f {config[filter_reads]}"
        " {input.dbs_fasta}"
        " {input.abc_fastas} > {log}"

rule copy_report:
    output:
        "report.ipynb"
    run:
        import pkg_resources
        report_path = pkg_resources.resource_filename("dbspro", 'report_template.ipynb')
        shell("jupyter nbconvert --to notebook {report_path} --output {output} --output-dir .")


rule make_report:
    output:
          "report.html"
    input:
         nb="report.ipynb",
         umis="umi_counts.tsv",
         reads="read_counts.tsv"
    shell:
         """
         jupyter nbconvert --execute --to notebook --inplace {input.nb}
         jupyter nbconvert --to html {input.nb}
         """

