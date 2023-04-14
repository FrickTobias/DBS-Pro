import pandas as pd
from snakemake.utils import validate

from dbspro.utils import get_abcs
from dbspro.cli.init import CONFIGURATION_FILE_NAME, ABC_FILE_NAME, SAMPLE_FILE_NAME

# Read sample and handles files.
configfile: CONFIGURATION_FILE_NAME
validate(config, "config.schema.yaml")

abc = get_abcs(ABC_FILE_NAME)
samples = pd.read_csv(SAMPLE_FILE_NAME, sep="\t").set_index("Sample", drop=False)

# Get required values
abc_len = len(abc["Sequence"][0]) - 1
dbs_n = "N"*len(config["dbs"])


rule all:
    input: 'report.html', 'data.tsv.gz'


if config["h1"] is None: # For PBA input
    dbs_trim = f"-a {config['h2']}"
    abs_umi_adapter = f"^{dbs_n}{config['h2']}...{config['h3']}"
else: # For DBS-Pro input
    dbs_trim = f"-g ^{config['h1']}...{config['h2']}"
    abs_umi_adapter = f"^{config['h1']}{dbs_n}{config['h2']}...{config['h3']}"

do_sampling = "subsampled." if config["subsample"] != -1 else ""
nr_samples = len(samples)


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


rule extract_dbs:
    """Extract DBS and trim handle between DBS and ABC."""
    output:
        reads="{sample}.dbs-raw.fastq.gz"
    input:
        reads=f"{{sample}}.{do_sampling}fastq.gz"
    log: "log_files/{sample}.cutadapt-extract-dbs.log"
    threads: max(workflow.cores / nr_samples, 4)
    params:
        trim=dbs_trim,
        err_rate=config["trim_err_rate"],
        min_len=len(config["dbs"]) - config["dbs_len_span"],
        max_len=len(config["dbs"]) + config["dbs_len_span"],
    shell:
        "cutadapt"
        " {params.trim}"
        " --discard-untrimmed"
        " -e {params.err_rate}"
        " -m {params.min_len}"
        " -M {params.max_len}"
        " --max-n 0"
        " -j {threads}"
        " -o {output.reads}"
        " {input.reads}"
        " > {log}"


rule extract_abc_umi:
    """Extract ABC and UMI."""
    output:
        reads="{sample}.trimmed-abc.fastq.gz"
    input:
        reads=f"{{sample}}.{do_sampling}fastq.gz"
    log: "log_files/{sample}.cutadapt-extract-abc-umi.log"
    threads: max(workflow.cores / nr_samples, 4)
    params:
        adapter=abs_umi_adapter,
        err_rate=config["trim_err_rate"],
        min_len=abc_len + config["umi_len"] - config["abc_umi_len_span"],
        max_len=abc_len + config["umi_len"] + config["abc_umi_len_span"],
    shell:
        "cutadapt"
        " -g {params.adapter}"
        " --discard-untrimmed"
        " -e {params.err_rate}"
        " -m {params.min_len}"
        " -M {params.max_len}"
        " -j {threads}"
        " -o {output.reads}"
        " --max-n 0"
        " {input.reads}"
        " > {log}"


rule demultiplex_abc:
    """Demultiplexes ABC sequnces and trims it to give ABC-specific UMI fastq files."""
    output:
        reads=touch(expand("ABCs/{{sample}}.{name}-UMI-raw.fastq.gz", name=abc['Target']))
    input:
        reads="{sample}.trimmed-abc.fastq.gz"
    log: "log_files/{sample}.cutadapt-id-abc.log"
    params:
        file=config["abc_file"],
        err_rate=config["demultiplex_err_rate"]
    shell:
        "cutadapt"
        " -g file:{params.file}"
        " --no-indels"
        " -e {params.err_rate}"
        " -o ABCs/{wildcards.sample}.{{name}}-UMI-raw.fastq.gz"
        " {input.reads}"
        " > {log}"


rule dbs_cluster:
    """Cluster DBS sequence using starcode for error correction."""
    output:
        clusters="{sample}.dbs-clusters.txt.gz"
    input:
        reads="{sample}.dbs-raw.fastq.gz"
    log: "log_files/{sample}.starcode-dbs-cluster.log"
    threads: max(workflow.cores / nr_samples, 4)
    params:
        dist = config["dbs_cluster_dist"]
    shell:
        "pigz -cd {input.reads} |"
        " starcode"
        " --print-clusters"
        " -t {threads}"
        " -d {params.dist}"
        " 2> {log} | pigz > {output.clusters}"


rule umi_cluster:
    """Cluster UMIs using UMI-tools API for each DBS and ABC to error correct them."""
    output:
        reads="ABCs/{sample}.{target}-UMI-corrected.fasta.gz"
    input:
        abc_reads="ABCs/{sample}.{target}-UMI-raw.fastq.gz",
        dbs_corrected="{sample}.dbs-corrected.fasta.gz"
    log: "log_files/{sample}.splitcluster-{target}.log"
    params:
        dist = config["abc_cluster_dist"],
        length = config["umi_len"]
    shell:
        "dbspro splitcluster"
        " {input.dbs_corrected}"
        " {input.abc_reads}"
        " -o {output.reads}"
        " -t {params.dist}"
        " -l {params.length}"
        " 2> {log}"


rule correct_dbs:
    """Combine DBS clustering results with original FASTA for error correction."""
    output:
        reads="{sample}.dbs-corrected.fasta.gz"
    input:
        reads="{sample}.dbs-raw.fastq.gz",
        clusters="{sample}.dbs-clusters.txt.gz"
    log: "log_files/{sample}.correctfastq-dbs.log"
    params:
        dbs = config['dbs']
    shell:
        "dbspro correctfastq"
        " {input.reads}"
        " {input.clusters}"
        " --output-fasta {output.reads}"
        " --barcode-pattern {params.dbs}"
        " 2> {log}"


rule integrate:
    """Integrate data into TSV with each DBS, ABC, UMI combination with read count for each sample."""
    output:
        data="{sample}.data.tsv.gz"
    input:
        dbs_fasta="{sample}.dbs-corrected.fasta.gz",
        abc_fastas=expand("ABCs/{{sample}}.{abc}-UMI-corrected.fasta.gz", abc=abc['Target'])
    log: "log_files/{sample}.integrate.log"
    shell:
        "dbspro integrate"
        " -o {output.data}"
        " {input.dbs_fasta}"
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


rule make_report:
    """Make jupyter notebook for final report"""
    output:
          html = "report.html",
          notebook = "report.ipynb",
    input:
         data="data.tsv.gz"
    log: "log_files/make_report.log"
    run:
        import pkg_resources
        report_path = pkg_resources.resource_filename("dbspro", 'report_template.ipynb')
        shell(
            "jupyter nbconvert --ClearMetadataPreprocessor.enabled=True --to notebook {report_path} --output {output.notebook} --output-dir . 2>> >(tee {log} >&2);"
            " jupyter nbconvert --execute --to notebook --inplace {output.notebook} 2>> >(tee {log} >&2);"
            " jupyter nbconvert --to html {output.notebook} 2>> >(tee {log} >&2)"
         )
