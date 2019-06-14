#!/bin/bash
set -euo pipefail


threads="nej!"


if  [[ "Linux" == $(uname) ]]
then
    threads=nproc
elif [[ "Darwin" == $(uname) ]]
then
    threads=$(sysctl -n hw.ncpu)
fi

# Argparsing
while getopts "t:h" OPTION
do
    case ${OPTION} in

        t)
            threads=${OPTARG}
            ;;
        h)
            printf 'DBSpro_automation.sh

Useage:     bash DBSpro_automation.sh <options> <reads.fq> <output_dir>
NB:         options must be given before arguments.

Positional arguments (REQUIRED)
  <reads.fq>    Read one in .fastq format. Also handles gzip files (.fastq.gz)
  <output_dir>  Output directory for analysis results

Global optional arguments
  -t  numbers of processors for threading                                               DEFAULT: 1
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
    exit 0e
fi

# PATH to WGH_Analysis folder
DBSpro_path=$(dirname "$0")

# output folder
path=$ARG2
mkdir -p $path

ln -s $PWD/$ARG1 $path/reads.fastq.gz
snakemake -j $threads $path/umi-counts.txt $path/umi-density-plot.png $path/read-density-plot.png