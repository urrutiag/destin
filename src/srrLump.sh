#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH --mem=32g

while getopts f:i:c:o: opt; do
  case ${opt} in
    f ) SRRFile=$OPTARG  ;;
    i ) start=$OPTARG ;;
    c ) count=$OPTARG ;;
    o ) outDir=$OPTARG
  esac
done

SRRNames=( $( cat $SRRFile ) )
echo "${SRRNames[@]:$start:$count}"

fastq-dump "${SRRNames[@]:$start:$count}" --split-files -O $outDir/fastq --gzip

