# Cutadapt trimming

ABC_list = ["ABC1-GCGTA" , "ABC2-ATAGC", "ABC3-GTGCA"]

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
        reads="{dir}/dbs.fastq.gz"
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

rule identify_abc_1:
    "Identifies ABC and trims it to give ABC-specific UMI fastq files."
    output:
        reads="{dir}/ABC1-GCGTA-UMI.fastq.gz"
    input:
        reads="{dir}/trimmed-abc.fastq.gz"
    log: "{dir}/id-abc1.log"
    threads: 20
    shell:
        "cutadapt"
        " -g ^GCGTA"
        " -m 6"
        " -M 6"
        " -e 0.2"
        " -j {threads}"
        " -o {output.reads}"
        " {input.reads}"
        " > {log}"

rule identify_abc_2:
    "Identifies ABC and trims it to give ABC-specific UMI fastq files."
    output:
        reads="{dir}/ABC2-ATAGC-UMI.fastq.gz"
    input:
        reads="{dir}/trimmed-abc.fastq.gz"
    log: "{dir}/id-abc2.log"
    threads: 20
    shell:
        "cutadapt"
        " -g ^ATAGC"
        " -m 6"
        " -M 6"
        " -e 0.2"
        " -j {threads}"
        " -o {output.reads}"
        " {input.reads}"
        " > {log}"

rule identify_abc_3:
    "Identifies ABC and trims it to give ABC-specific UMI fastq files."
    output:
        reads="{dir}/ABC3-GTGCA-UMI.fastq.gz"
    input:
        reads="{dir}/trimmed-abc.fastq.gz"
    log: "{dir}/id-abc3.log"
    threads: 20
    shell:
        "cutadapt"
        " -g ^GTGCA"
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
        reads="{dir}/dbs.fastq.gz"
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
        clusters=expand("{{dir}}/{ABC}-UMI-clusters.txt", ABC=ABC_list)
    input:
        reads=expand("{{dir}}/{ABC}-UMI.fastq.gz", ABC=["ABC1-GCGTA" , "ABC2-ATAGC", "ABC3-GTGCA"])
    log: "{dir}/abc-clusters.log"
    threads: 20
    shell:
        "pigz -cd {input.reads} | starcode"
        " --print-clusters"
        " -t {threads}"
        " -d 1"
        " -o {output.clusters}"

# DBSpro

rule error_correct:
    "Combine cluster results with original files to error correct them."
    output:
        reads=expand("{{dir}}/{corr_file}-corrected.fastq.gz", corr_file=["dbs", "ABC1-GCGTA-UMI", "ABC2-ATAGC-UMI", "ABC3-GTGCA-UMI"])
    input:
        reads=expand("{{dir}}/{corr_file}.fastq.gz", corr_file=["dbs", "ABC1-GCGTA-UMI", "ABC2-ATAGC-UMI", "ABC3-GTGCA-UMI"]),
        clusters=expand("{{dir}}/{corr_file}-clusters.txt", corr_file=["dbs", "ABC1-GCGTA-UMI", "ABC2-ATAGC-UMI", "ABC3-GTGCA-UMI"])
    log: "{dir}/error-correct.log"
    threads: 20
    shell:
        "DBSpro correctfastq"
        " input.reads"
        " input.clusters"
        " output.reads"

rule analyze:
    "Analyzes all result files"
    output:
        counts="{dir}/umi-counts.txt",
        umi_plot="{dir}/umi-density-plot.png",
        reads_plot="{dir}/read-density-plot.png"
    input:
        dbs_fastq="{dir}/dbs-corrected.fastq.gz",
        abc1_fastq="{dir}/ABC1-GCGTA-UMI-corrected.fastq.gz",
        abc2_fastq="{dir}/ABC2-ATAGC-UMI-corrected.fastq.gz",
        abc3_fastq="{dir}/ABC3-GTGCA-UMI-corrected.fastq.gz"
    log: "{dir}/analyze.log"
    threads: 20
    shell:
        "DBSpro analyze"
        " -f 4"
        " input.dbs_fastq"
        " input.abc1_fastq"
        " input.abc2_fastq"
        " input.abc3_fastq"
        " output.counts"
        " output.umi_plot"
        " output.reads.plot"