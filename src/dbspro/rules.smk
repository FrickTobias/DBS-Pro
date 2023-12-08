"""
Snakefile for DBS-Pro pipeline
"""
from importlib.resources import files, as_file

import pandas as pd
from snakemake.utils import validate

from dbspro.utils import get_abcs
from dbspro.cli.init import CONFIGURATION_FILE_NAME, ABC_FILE_NAME, SAMPLE_FILE_NAME, MULTIQC_CONFIG_NAME

# Read sample and handles files.
configfile: CONFIGURATION_FILE_NAME
validate(config, "config.schema.yaml")

abc = get_abcs(ABC_FILE_NAME)
samples = pd.read_csv(SAMPLE_FILE_NAME, sep="\t").set_index("Sample", drop=False)

# Get required values
abc_len = len(abc["Sequence"].iloc[0]) - 1
abc_umi_len = abc_len + config["umi_len"]
dbs_n = "N"*len(config["dbs"])
dbs_h2_abs_umi_len = len(config["dbs"]) + len(config["h2"]) + abc_len + config["umi_len"]

if config["h1"] is None: # For PBA input
    trim_outer = f"-a {config['h3']}"
    abs_umi_adapter = f"^{dbs_n}{config['h2']}...{config['h3']}"
else: # For DBS-Pro input
    trim_outer = f"-g ^{config['h1']}...{config['h3']}"
    abs_umi_adapter = f"^{config['h1']}{dbs_n}{config['h2']}...{config['h3']}"

do_sampling = "subsampled." if config["subsample"] != -1 else ""
nr_samples = len(samples)


wildcard_constraints:
    sample="\w+"


rule all:
    input: 
        'report.html', 
        'data.tsv.gz', 
        'multiqc_report.html',
        expand("{sample}.counts.h5ad", sample=samples["Sample"])


rule subsample:
    """Subsample input reads to required depth if needed"""
    output:
        reads = "{sample}.subsampled.fastq.gz"
    input:
        reads = "{sample}.fastq.gz"
    params:
        number = config["subsample"] if config["subsample"] > 0 else samples["Reads"].min()
    threads: max(workflow.cores / nr_samples, 4)
    run:
        # Only subsample file if needed
        if samples.loc[wildcards.sample, "Reads"] > params.number:
            shell(
                "seqtk sample"
                " -s 9999"
                " {input.reads}"
                " {params.number}"
                " | pigz > {output.reads}"
            )
        else:
            shell("ln -s {input.reads} {output.reads}")


rule fastqc:
    output:
        html="log_files/{sample}_fastqc.html",
        zip="log_files/{sample}_fastqc.zip"
    input:
        reads=f"{{sample}}.{do_sampling}fastq.gz"
    log: "log_files/{sample}_fastqc.log"
    shell:
        "fastqc"
        " -o log_files"
        " -t 1"
        " {input.reads}"
        " &> {log}"


rule trim_outer_handles:
    """Trim outer handles leaving DBS - H2 - ABC+UMI."""
    output:
        reads="{sample}.trimmed.fastq.gz",
    input:
        reads=f"{{sample}}.{do_sampling}fastq.gz"
    log: "log_files/{sample}.trimmed.log"
    threads: max(workflow.cores / nr_samples, 4)
    params:
        trim=trim_outer,
        err_rate=config["trim_err_rate"],
        min_len=dbs_h2_abs_umi_len - int(dbs_h2_abs_umi_len * 0.1),
        max_len=dbs_h2_abs_umi_len + int(dbs_h2_abs_umi_len * 0.1),
        overlap=5
    shell:
        "cutadapt"
        " {params.trim}"
        " -e {params.err_rate}"
        " -m {params.min_len}"
        " -M {params.max_len}"
        " --max-n 0"
        " -O {params.overlap}"
        " -Z"
        " -j {threads}"
        " -o {output.reads}"
        " {input.reads}"
        " > {log}"


rule fastqc_trimmed:
    output:
        html="log_files/{sample}.trimmed_fastqc.html",
        zip="log_files/{sample}.trimmed_fastqc.zip"
    input:
        reads="{sample}.trimmed.fastq.gz"
    log: "log_files/{sample}.trimmed_fastqc.log"
    shell:
        "fastqc"
        " -o log_files"
        " -t 1"
        " {input.reads}"
        " &> {log}"


rule extract_dbs_abc_umi:
    """Extract DBS and ABC+UMI."""
    output:
        dbs="{sample}.trimmed.dbs.fasta.gz",
        abc_umi="{sample}.trimmed.abc_umi.fasta.gz",
        dbs_tmp=temp("{sample}.trimmed.dbs.txt"),
    input:
        reads="{sample}.trimmed.fastq.gz"
    log: "log_files/{sample}.trimmed.dbs.log",
    threads: max(workflow.cores // nr_samples, 4)
    params:
        err_rate=config["trim_err_rate"],
        dbs_len = len(config["dbs"]),
        abs_umi_len = abc_len + config["umi_len"],
        handle = dbs_n + config["h2"],
    shell:
        # Extract DBS and ABC+UMI using cutadapt
        "cutadapt"
        " -g ^{params.handle}"
        " -j {threads}"
        " --rename {{id}}"
        " -m {params.abs_umi_len}"
        " -M {params.abs_umi_len}"
        " --fasta"
        " --discard-untrimmed"
        " -o {output.abc_umi}"
        " --wildcard-file {output.dbs_tmp}"
        " -Z"
        " {input.reads}"
        " > {log}"
        " && "
        # Convert TXT file with DBS sequences to FASTA
        "awk -F' ' '{{print \">\"$2\"\\n\"$1 }}' < {output.dbs_tmp}"
        " | "
        "pigz -1 > {output.dbs}"


rule count_dbs:
    """Count DBS sequences. This saves memory for the starcode step."""
    input:
        reads="{sample}.trimmed.dbs.fasta.gz"
    output:
        counts=temp("{sample}.trimmed.dbs.counts.tsv")
    shell:
        "gunzip -c {input.reads}"
        " | "
        "awk -v OFS=\"\t\" '{{ if ( NR%2==0 ) {{ counts[$1]++ }} }} END {{ for (barcode in counts) print (barcode,counts[barcode]) }}'"
        " > {output.counts}"


rule dbs_cluster:
    """Cluster DBS sequence using starcode for error correction."""
    output:
        clusters="{sample}.trimmed.dbs.clusters.txt.gz"
    input:
        reads="{sample}.trimmed.dbs.counts.tsv"
    log: "log_files/{sample}.dbs.clusters.log"
    threads: max(workflow.cores / nr_samples, 4)
    params:
        dist = config["dbs_cluster_dist"]
    shell:
        "starcode"
        " --print-clusters"
        " -t {threads}"
        " -d {params.dist}"
        " {input.reads}"
        " 2> {log} | pigz > {output.clusters}"


rule correct_dbs:
    """Combine DBS clustering results with original FASTA for error correction."""
    output:
        reads="{sample}.trimmed.dbs.corrected.fasta.gz"
    input:
        reads="{sample}.trimmed.dbs.fasta.gz",
        clusters="{sample}.trimmed.dbs.clusters.txt.gz"
    log: "log_files/{sample}.trimmed.dbs.corrected.log"
    shell:
        "dbspro correctfastq"
        " {input.reads}"
        " {input.clusters}"
        " --output-fasta {output.reads}"
        " 2> {log}"


rule tagfastq:
    """Tag ABC and UMI sequences with DBS sequence and sort by barcode."""
    output:
        reads="{sample}.trimmed.abc_umi.tagged.fasta.gz"
    input:
        dbs="{sample}.trimmed.dbs.corrected.fasta.gz",
        abc_umi="{sample}.trimmed.abc_umi.fasta.gz"
    log: "log_files/{sample}.trimmed.abc_umi.tagged.log"
    shell:
        "dbspro tagfastq"
        " {input.abc_umi}"
        " {input.dbs}"
        " -s ' '"
        " 2> {log}"
        " | "
        "paste - -"
        " | "
        "awk -F ' ' '{{OFS=\"\t\"; print $2,$0}}'"
        " | "
        "sort -t $'\t' -k1,1"
        " | "
        "cut -f 2-"
        " | "
        "tr '\t' '\n'"
        " | "
        "pigz -1 > {output.reads}"


rule demultiplex_abc:
    """Demultiplexes ABC sequnces and trims it to give ABC-specific UMI fastq files."""
    output:
        reads=touch(expand("ABCs/{{sample}}.{name}.umi.fasta.gz", name=abc['Target']))
    input:
        reads="{sample}.trimmed.abc_umi.tagged.fasta.gz"
    log: 
        log = "log_files/{sample}.abc.umi.log",
        json = "log_files/{sample}.abc.umi.json"
    params:
        file=config["abc_file"],
        err_rate=config["demultiplex_err_rate"]
    shell:
        "cutadapt"
        " -g file:{params.file}"
        " --no-indels"
        " -e {params.err_rate}"
        " --json {log.json}"
        " -Z"
        " -o ABCs/{wildcards.sample}.{{name}}.umi.fasta.gz"
        " {input.reads}"
        " > {log.log}"


rule umi_cluster:
    """Cluster UMIs using UMI-tools API for each DBS and ABC to error correct them."""
    output:
        reads="ABCs/{sample}.{target}.umi.corrected.fasta.gz"
    input:
        reads="ABCs/{sample}.{target}.umi.fasta.gz",
    log: "log_files/{sample}.{target}.umi.corrected.log"
    params:
        dist = config["abc_cluster_dist"],
        length = config["umi_len"]
    run:
        if params.dist > 0:
            shell(
                "dbspro splitcluster"
                " {input.reads}"
                " -o {output.reads}"
                " -t {params.dist}"
                " -l {params.length}"
                " 2> {log}"
            )
        else: # Copy file if no clustering specified
            shell(
                "cp {input.reads} {output.reads}"
            )


rule integrate:
    """Integrate data into TSV with each DBS, ABC, UMI combination with read count for each sample."""
    output:
        data="{sample}.data.tsv.gz"
    input:
        abc_fastas=expand("ABCs/{{sample}}.{abc}.umi.corrected.fasta.gz", abc=abc['Target'])
    log: "log_files/{sample}.integrate.log"
    params:
        dbs = config['dbs']
    shell:
        "dbspro integrate"
        " -o {output.data}"
        " --barcode-pattern {params.dbs}"
        " {input.abc_fastas}"
        " 2> {log}"


rule merge_data:
    """Merge data from all samples"""
    output:
        data = "data.tsv.gz"
    input: 
        data_files = expand("{sample}.data.tsv.gz", sample=samples["Sample"])
    run:
        merged_data = pd.concat([pd.read_csv(file, sep="\t") for file in input.data_files])
        merged_data.to_csv(output.data, sep="\t", index=False)


rule generate_h5ad:
    """Generate h5ad file for downstream analysis using scanpy"""
    output:
        data = "{sample}.counts.h5ad"
    input:
        data = "{sample}.data.tsv.gz"
    threads: workflow.cores
    script: "scripts/generate_h5ad.py"


rule get_preseq_vals:
    """Prepare for running preseq"""
    input:
        data = "{sample}.data.tsv.gz"
    output:
        txt = temp("{sample}.vals.txt")
    shell:
        "zless {input.data} | tail -n+2 | cut -f 4 > {output.txt}"


rule preseq:
    """Run preseq to estimate library complexity"""
    input:
        txt = "{sample}.vals.txt"
    output:
        txt = "log_files/{sample}.preseq.txt"
    shell:
        "preseq lc_extrap"
        " -o {output.txt}"
        " -V"
        " -e 1e+07"
        " {input.txt}"

rule preseq_real_counts:
    """TSV with real counts for preseq MultiQC plot"""
    input:
        data = expand("{sample}.data.tsv.gz", sample=samples["Sample"])
    output:
        tsv = "preseq_real_counts.tsv"
    run:
        with open(output.tsv, "w") as f:
            # Columns are: Sample name, number of reads, number of unique constructs.
            for file in input.data:
                sample = file.split(".")[0]
                df = pd.read_csv(file, sep="\t")
                reads = df["ReadCount"].sum()
                constructs = len(df)
                f.write(f"{sample}\t{reads}\t{constructs}\n")


rule make_report:
    """Make jupyter notebook for final report"""
    output:
          html = "report.html",
          notebook = "report.ipynb",
    input:
         data="data.tsv.gz"
    log: "log_files/make_report.log"
    run:
        with as_file(files("dbspro").joinpath("report_template.ipynb")) as report_path:
            shell(
                "jupyter nbconvert --ClearMetadataPreprocessor.enabled=True --to notebook {report_path} --output {output.notebook} --output-dir . 2>> >(tee {log} >&2);"
                " jupyter nbconvert --execute --to notebook --inplace {output.notebook} 2>> >(tee {log} >&2);"
                " jupyter nbconvert --to html {output.notebook} 2>> >(tee {log} >&2)"
            )


def get_multiqc_config():
    """Get config for multiqc"""
    # Use user config if exists in workdir, otherwise use default config
    if os.path.exists(MULTIQC_CONFIG_NAME):
        return MULTIQC_CONFIG_NAME
    else:
        return files("dbspro").joinpath(MULTIQC_CONFIG_NAME)


rule multiqc:
    """Make multiqc report"""
    output:
        html="multiqc_report.html",
        dir=directory("multiqc_data")
    input:
        expand(rules.fastqc.output.zip, sample=samples["Sample"]),
        expand(rules.fastqc_trimmed.output.zip, sample=samples["Sample"]),
        expand(rules.trim_outer_handles.output.reads, sample=samples["Sample"]),
        expand(rules.preseq.output.txt, sample=samples["Sample"]),
        rules.preseq_real_counts.output.tsv,
    params:
        config=get_multiqc_config(),
    shell:
        "multiqc -c {params.config} -f log_files &> multiqc_report.html.log"
