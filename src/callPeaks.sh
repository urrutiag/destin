#!/bin/bash

#SBATCH --mem=32g
#SBATCH --time=24:00:00

while getopts s:o:l:q: opt; do
  case ${opt} in
    s ) sampleName=$OPTARG  ;;
    o ) outDir=$OPTARG ;;
    l ) blacklistFile=$OPTARG ;;
    q ) genomeMacs2=$OPTARG
  esac
done

bamDir=$outDir/bam
peaksDir=$outDir/peaks

finalBams=($bamDir/**.bam)

echo sampleName: $sampleName
echo PeaksDir: $peaksDir
echo bamDir: $bamDir
echo blacklistFile: $blacklistFile
echo genomeMacs2: $genomeMacs2
echo first10Bams: ${finalBams[@]:0:10}

#MACS2
#-g hs or mm#remove blacklist
macs2 callpeak -t "${finalBams[@]}" \
  -f BAM -g $genomeMacs2 \
  -n $peaksDir/"$sampleName"   \
  --nomodel  -p 0.01

#remove blacklist
bedtools subtract -a $peaksDir/"$sampleName"_peaks.narrowPeak \
  -b $blacklistFile > \
  $peaksDir/"$sampleName"_peaks.blacklist_removed.narrowPeak

