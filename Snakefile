# Cutadapt trimming


rule trim_3prime:
    "Trim 3' end for coupling sequence TTATATCACGACAAGAG."
    output:
        reads="{dir}/trimmed-3prim.fastq.gz"
    input:
        reads="{dir}/reads.fastq.gz"
    log: "{dir}/reads-3prim.log"
    threads: 20
    shell:
        "cutadapt"
        " -a TTATATCACGACAAGAG"
        " --discard-untrimmed"
        " -e 0.2"
        " -j {threads}"
        " -o {output.reads}"
        " {input.reads}"
        " > {log}"

rule extract_dbs:
    "Extract DBS and trim handle between DBS and ABC. H1+H2: CGATGCTAATCAGATCA, H3: AAGAGTCAATAGACCAT, H4: CTAACAGGATTCAGGTA"
    output:
        reads="{dir}/dbs-raw.fastq.gz"
    input:
        reads="{dir}/trimmed-3prim.fastq.gz"
    log: "{dir}/extract-dbs.log"
    threads: 20
    shell:
        "cutadapt"
        " -g ^CGATGCTAATCAGATCA...AAGAGTCAATAGACCATCTAACAGGATTCAGGTA"
        " --discard-untrimmed"
        " -e 0.2"
        " -j {threads}"
        " -o {output.reads}"
        " {input.reads}"
        " > {log}"

rule trim_to_abc:
    "Trim 3' end for coupling sequence between DBS and ABC."
    output:
        reads="{dir}/trimmed-abc.fastq.gz"
    input:
        reads="{dir}/trimmed-3prim.fastq.gz"
    log: "{dir}/reads-3prim.log"
    threads: 20
    shell:
        "cutadapt"
        " -g AAGAGTCAATAGACCATCTAACAGGATTCAGGTA"
        " --discard-untrimmed"
        " -e 0.2"
        " -j {threads}"
        " -o {output.reads}"
        " {input.reads}"
        " > {log}"

## ABCs

rule identify_abc:
    "Identifies ABC and trims it to give ABC-specific UMI fastq files."
    output:
        reads="{dir}/{abc}-UMI-raw.fastq.gz"
    input:
        reads="{dir}/trimmed-abc.fastq.gz"
    log: "{dir}/id-{abc}.log"
    threads: 20
    shell:
        "cutadapt"
        " -g ^{wildcards.abc}"
        " -m 6"
        " -M 6"
        " -e 0.2"
        " -j {threads}"
        " -o {output.reads}"
        " {input.reads}"
        " > {log}"

# Starcode clustering

rule dbs_cluster:
    "Cluster DBS sequence using starcode"
    output:
        clusters="{dir}/dbs-clusters.txt"
    input:
        reads="{dir}/dbs-raw.fastq.gz"
    log: "{dir}/dbs-clusters.log"
    threads: 20
    shell:
        "pigz -cd {input.reads} | starcode"
        " --print-clusters"
        " -t {threads}"
        " -d 2"
        " -o {output.clusters}"

rule abc_cluster:
    "Cluster ABC sequence using starcode"
    output:
        "{dir}/{sample}-UMI-clusters.txt"
    input:
        "{dir}/{sample}-UMI-raw.fastq.gz"
    log: "{dir}/{sample}-clusters.log"
    threads: 20
    shell:
        "pigz -cd {input} | starcode"
        " --print-clusters"
        " -t {threads}"
        " -d 1"
        " -o {output}"

# DBSpro

rule error_correct:
    "Combine cluster results with original files to error correct them."
    output:
        reads="{dir}/{corr_file}-corrected.fastq"
    input:
        reads="{dir}/{corr_file}-raw.fastq.gz",
        clusters="{dir}/{corr_file}-clusters.txt"
    log: "{dir}/{corr_file}error-correct.log"
    threads: 20
    shell:
        "DBSpro correctfastq"
        " {input.reads}"
        " {input.clusters}"
        " {output.reads}"

rule analyze:
    "Analyzes all result files"
    output:
        counts="{dir}/umi-counts.txt",
        umi_plot="{dir}/umi-density-plot.png",
        reads_plot="{dir}/read-density-plot.png"
    input:
        dbs_fastq="{dir}/dbs-corrected.fastq",
        abc1_fastq="{dir}/GCGTA-UMI-corrected.fastq",
        abc2_fastq="{dir}/ATAGC-UMI-corrected.fastq",
        abc3_fastq="{dir}/GTGCA-UMI-corrected.fastq"
    log: "{dir}/analyze.log"
    threads: 20
    shell:
        "DBSpro analyze"
        " -f 4"
        " {input.dbs_fastq}"
        " {input.abc1_fastq}"
        " {input.abc2_fastq}"
        " {input.abc3_fastq}"
        " {output.counts}"
        " {output.umi_plot}"
        " {output.reads_plot}"