#!/bin/bash

#SBATCH --mem=32g
#SBATCH --time=24:00:00

while getopts s:o:g:l:q: opt; do
  case ${opt} in
    s ) sampleName=$OPTARG  ;;
    o ) outDir=$OPTARG ;;
    g ) genomeBedtools=$OPTARG ;;
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
echo genomeBedtools: $genomeBedtools
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

#set start and end to summit
# awk 'BEGIN{ OFS="\t" } { summit=$10;  print $1, $2+$10, $2+$10, $4, $5, $6, $7, $8, $9, $10 }' \
#  $peaksDir/"$sampleName"_peaks.blacklist_removed.narrowPeak > \
#  $peaksDir/"$sampleName".summit.narrowPeak

#expand to +/- 250 from summit
#bedtools slop -i $peaksDir/"$sampleName".summit.narrowPeak \
#  -g $genomeBedtools \
#  -b 250 > $peaksDir/"$sampleName".500.narrowPeak

