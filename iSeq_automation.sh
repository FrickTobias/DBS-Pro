#! /bin/bash
#set -euo pipefail

processors=1
mailing=false

# Argparsing
while getopts "m:p:h" OPTION
do
    case ${OPTION} in

        p)
            processors=${OPTARG}
            ;;
        m)
            email=${OPTARG}
            mailing=true
            ;;
        h)
            printf 'iSeq_automation.sh

Useage:     bash iSeq_automation.sh <options> <reads.fq> <output_dir>
NB:         options must be given before arguments.

Pipeline outline:
  0.
  1.
  2.
  3.

Positional arguments (REQUIRED)
  <reads.fq>    Read one in .fastq format. Also handles gzip files (.fastq.gz)
  <output_dir>  Output directory for analysis results

Global optional arguments
  -m  mails the supplied email when analysis is finished                                DEFAULT: None
  -p  numbers of processors for threading                                               DEFAULT: 1
  -h  help (this output)                                                                DEFAULT: N/A
  \n'
	        exit 0
	        ;;
    esac
done

# Positonal redundancy for option useage
ARG1=${@:$OPTIND:1}
ARG2=${@:$OPTIND+1:1}

# Error handling
if [ -z "$ARG1" ] || [ -z "$ARG2" ]
then
    echo ""
    echo "ARGUMENT ERROR"
    echo "Did not find all positional arguments, see -h for more information."
    echo "(got reads:"$ARG1" and output:"$ARG2" instead)"
    echo ""
    exit 0
fi

printf '\n0. Argparsing & options'
printf '\nReads:\t'$ARG1'\nOutput:\t'$ARG2'\n'
printf '\nThreads:\t'$processors

# Mailing option
if $mailing
then
    if [[ $email == *"@"* ]]
    then
        printf '\nMail:\t\t'$email
    else
        echo ''
        echo 'OPTION ERROR: -m '
        echo ''
        echo 'Please supply email on format john.doe@domain.org'
        echo '(got "'$email'" instead)'
        exit 0
    fi
fi

# PATH to WGH_Analysis folder
iSeq_path=$(dirname "$0")

# output folder
path=$ARG2
mkdir -p $path

# File one prep
file=$ARG1
name_ext=$(basename "$file")
name="${name_ext%.*}"
file_name="$path/${name_ext%.*}"

# Logfiles
trim_logfile=$path'/1_trim.log'
bc_cluster_logfile=$path'/2_cluster.log'
UMI_cluster_logfile=$path'/3_map.log'
abc_identification_logfile=$path'/4_rmdup.log'

# Mailing
if $mailing
then
    echo 'ANALYSIS STARTING '$(date) | mail -s $path $email
fi

printf '\n\n'"`date`"'\tANALYSIS STARTING\n'


#
# # # 1. Read trimming
#


# 5' trimming for H5
# H5: TTATATCACGACAAGAG
cutadapt -a TTATATCACGACAAGAG \
    --discard-untrimmed \
    -e 0.2 \
    -j $processors \
    -o $file_name".5prim.fastq.gz" \
    $ARG1 > $path/"5prim_trim.log"


# Extracting barcode located between H1-H2 and H3-H4 (H1-H2-DBS-H3-H4)
# H1+H2: CGATGCTAATCAGATCA
# H3: AAGAGTCAATAGACCAT
# H4: CTAACAGGATTCAGGTA
cutadapt -g ^CGATGCTAATCAGATCA...AAGAGTCAATAGACCATCTAACAGGATTCAGGTA \
    --discard-untrimmed \

    -e 0.2 \
    -j $processors \
    -o $file_name".DBS.fastq.gz" \
    $file_name".5prim.fastq.gz"  > $path/"bc_trim.log"

# Trimming ABC/UMI file
# Extracting barcode located between H1-H2 and H3-H4 (H1-H2-DBS-H3-H4-ABC-UMI-H5-i7)
cutadapt -g AAGAGTCAATAGACCATCTAACAGGATTCAGGTA \
    --discard-untrimmed \
    -e 0.2 \
    -j $processors \
    -o $file_name".ABC_UMI.fastq.gz" \
    $file_name".5prim.fastq.gz"  > $path/"ABC_UMI_trim.log"

# Identifying ABC using cutadapt to look for ABC:s
for abc in "GCGTA" "ATAGC" "GTGCA"
do
    cutadapt -g $abc \
        --discard-untrimmed \
        -e 0.2 \
        -m 6 \
        -M 6 \
        -j $processors \
        -o $file_name"."$abc".fastq.gz" \
        $file_name".ABC_UMI.fastq.gz" > $path/"ABC_trim."$abc".log"
done

# 2. Barcode and UMI error correction

# Barcode error correction
pigz -cd $file_name".DBS.fastq.gz" | starcode \
    --print-clusters \
    -t $processors \
    -d 4 \
    -o $file_name".DBS.SC_out"

# UMI error correction
for abc in "GCGTA" "ATAGC" "GTGCA"
do

    # If file is empty, skip file (= No UMI:s found for that ABC)
    ls -l $file_name"."$abc".fastq.gz"
    [ -s $file_name"."$abc".fastq.gz" ]
    empty=$(echo $?)

    if ! $empty
    then
        pigz -cd $file_name"."$abc".fastq.gz" | starcode \
            --print-clusters \
            -t $processors \
            -d 1 \
            -o $file_name"."$abc".SC_out"
    fi
done

# Incorporating the starcode clustering information in the fastq files
for err_corr_file in "DBS" "GCGTA" "ATAGC" "GTGCA"
do

    # If file is empty, skip file (= No UMI:s found for that ABC)
    ls -l $file_name"."$err_corr_file".fastq.gz"
    [ -s $file_name"."$err_corr_file".fastq.gz" ]
    empty=$(echo $?)

    if ! $empty
    then
        python3 $iSeq_path/python\ scripts/sum_analysis_prep.py \
            $file_name"."$err_corr_file".fastq.gz" \
            $file_name"."$err_corr_file".SC_out" \
            $file_name"."$err_corr_file".err_corr.fastq"
    fi
done

python3 $iSeq_path/python\ scripts/sum_results.py \
    $file_name".DBS.err_corr.fastq" \
    $file_name".GCGTA.err_corr.fastq" \
    $file_name".ATAGC.err_corr.fastq" \
    $file_name".GTGCA.err_corr.fastq" \
    $path"/result_summary.tsv"

























