import pandas as pd
import os

print(os.getcwd())
# Read sample and handles files.
abc = pd.read_csv(config["abc_sequences"], sep='\t').set_index("Antibody-target", drop=False)
handles = pd.read_csv(config["handles"], sep='\t').set_index("Name", drop=False)

# Get required values
abc_len = list(map(len, abc['Barcode-sequence']))[0]    # Assumes that all ABC are same length
abc_umi_len = abc_len + config["umi_len"]
dbs = "N"*config["dbs_len"]

rule all:
    input: 'umi-counts.txt', 'umi-density-plot.png', 'read-density-plot.png'

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
        " -g ^{handles[Sequence][h1]}...{handles[Sequence][h2]}"
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
        " -g ^{handles[Sequence][h1]}{dbs}{handles[Sequence][h2]}...{handles[Sequence][h3]}"
        " --discard-untrimmed"
        " -e 0.2"
        " -m {abc_umi_len}"
        " -M {abc_umi_len}"
        " -j {threads}"
        " -o {output.reads}"
        " {input.reads}"
        " > {log}"

## ABCs

"Identifies ABC and trims it to give ABC-specific UMI fastq files."
rule identify_abc:
    output:
        reads="{sample}-UMI-raw.fastq.gz"
    input:
        reads="trimmed-abc.fastq.gz"
    log: "log_files/cutadapt-id-abc-{sample}.log"
    threads: 20
    params:
        seq = lambda wildcards: abc['Barcode-sequence'][wildcards.sample]
    shell:
        "cutadapt"
        " -g ^{params.seq}"
        " -m {config[umi_len]}"
        " -M {config[umi_len]}"
        " -e 0.2"
        " -j {threads}"
        " -o {output.reads}"
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
        " -o {output.clusters}"


"Cluster ABC sequence using starcode"
rule abc_cluster:
    output:
        clusters="{sample}-UMI-clusters.txt"
    input:
        reads="{sample}-UMI-raw.fastq.gz"
    log: "log_files/starcode-abc-cluster-{sample}.log"
    threads: 20
    shell:
        "pigz -cd {input.reads} | starcode"
        " --print-clusters"
        " -t {threads}"
        " -d {config[abc_cluster_dist]}"
        " -o {output.clusters}"


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
        counts="umi-counts.txt",
        umi_plot="umi-density-plot.png",
        reads_plot="read-density-plot.png"
    input:
        dbs_fasta="dbs-corrected.fasta",
        abc_fastas=expand("{abc}-UMI-corrected.fasta", abc=abc['Antibody-target'])
    log: "log_files/analyze.log"
    threads: 20
    shell:
        "dbspro analyze"
        " -f {config[filter_reads]}"
        " {input.dbs_fasta}"
        " {output.counts}"
        " {output.umi_plot}"
        " {output.reads_plot}"
        " {input.abc_fastas}"
