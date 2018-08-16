#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH --mem=32g

while getopts f:i:c:o:s:b:r: opt; do
  case ${opt} in
    f ) cellIDFile=$OPTARG  ;;
    i ) start=$OPTARG ;;
    c ) count=$OPTARG ;;
    o ) outDir=$OPTARG ;;
    s ) sampleName=$OPTARG ;;
    b ) genomeBowtie=$OPTARG ;;
    r ) srcDir=$OPTARG
  esac
done

echo cellIDFile: $cellIDFile
echo start: $start
echo count: $count
echo outDir: $outDir
echo sampleName: $sampleName
echo genomeBowtie: $genomeBowtie
echo srcDir: $srcDir

cellIDs=( $( cat $cellIDFile ) )

#echo "${cellIDs[@]:$start:$count}"

 
for cellID in "${cellIDs[@]:$start:$count}"
do
# echo $cellID
  $srcDir/fastqToBam.sh \
    -c $cellID  \
    -s $sampleName \
    -o $outDir \
    -g $genomeBowtie
done

