rule trim_3prime:
    "Trim 3' end for coupling sequence TTATATCACGACAAGAG."
    output:
        reads=temp("{dir}/trimmed-3prim.fastq.gz")
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

#rule extract_dbs:
#    "Trim handles around dbs to identify sequence"
#    output:
#        reads=

#
## 5' trimming for H5
## H5: TTATATCACGACAAGAG
#
#
## Extracting barcode located between H1-H2 and H3-H4 (H1-H2-DBS-H3-H4)
## H1+H2: CGATGCTAATCAGATCA
## H3: AAGAGTCAATAGACCAT
## H4: CTAACAGGATTCAGGTA
#cutadapt -g ^CGATGCTAATCAGATCA...AAGAGTCAATAGACCATCTAACAGGATTCAGGTA \
#    --discard-untrimmed \
#    -e 0.2 \
#    -j $processors \
#    -o $file_name".DBS.fastq.gz" \
#    $file_name".5prim.fastq.gz"  > $path/"bc_trim.log"
#
## Trimming ABC/UMI file
## Extracting barcode located between H1-H2 and H3-H4 (H1-H2-DBS-H3-H4-ABC-UMI-H5-i7)
#cutadapt -g AAGAGTCAATAGACCATCTAACAGGATTCAGGTA \
#    --discard-untrimmed \
#    -e 0.2 \
#    -j $processors \
#    -o $file_name".ABC_UMI.fastq.gz" \
#    $file_name".5prim.fastq.gz"  > $path/"ABC_UMI_trim.log"
#
## Identifying ABC using cutadapt to look for ABC:s
#for abc in "GCGTA" "ATAGC" "GTGCA"
#do
#    cutadapt -g $abc \
#        --discard-untrimmed \
#        -e 0.2 \
#        -m 6 \
#        -M 6 \
#        -j $processors \
#        -o $file_name"."$abc".fastq.gz" \
#        $file_name".ABC_UMI.fastq.gz" > $path/"ABC_trim."$abc".log"
#done
#