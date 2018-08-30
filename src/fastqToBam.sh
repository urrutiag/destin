#!/bin/bash

barcodeSampleName="$barcode"."$sampleName"
while getopts c:s:o:g: opt; do
  case ${opt} in
    c ) cellID=$OPTARG  ;;
    s ) sampleName=$OPTARG  ;;
    o ) outDir=$OPTARG ;;
    g ) genomeBowtie=$OPTARG ;;
    w ) adaptor1=$OPTARG ;;
    x ) adaptor2=$OPTARG ;;
    y ) adaptor3=$OPTARG ;;
    z ) adaptor4=$OPTARG 
  esac
done


tempDir=$outDir/temp
bamDir=$outDir/bam
fastqDir=$outDir/fastq

echo cellID: $cellID
echo sampleName: $sampleName
echo fastqDir: $fastqDir
echo bamDir: $bamDir
echo tempDir: $tempDir
echo picardPath: $picardPath
echo genomeBowtie: $genomeBowtie

barcodeSampleName="$cellID"."$sampleName"

myFastqs=( $(ls $fastqDir | grep $cellID) )

#trim adaptor
cutadapt --minimum-length=20  \
  -A $adaptor1 -a $adaptor2 -G $adaptor3 \
  -g $adaptor4 -o $tempDir/$barcodeSampleName.R1.trimmed.fastq.gz -p \
  $tempDir/$barcodeSampleName.R2.trimmed.fastq.gz \
  $fastqDir/${myFastqs[0]} \
  $fastqDir/${myFastqs[1]}


#align
bowtie2  -x $genomeBowtie \
  -p 25 -t -X2000–no-mixed–no-discordant \
  -1 $tempDir/$barcodeSampleName.R1.trimmed.fastq.gz \
  -2 $tempDir/$barcodeSampleName.R2.trimmed.fastq.gz | samtools view -Sb /dev/stdin \
   > $tempDir/$barcodeSampleName.bam


#sort
java -Xmx4g -jar $picardPath/picard.jar SortSam \
  INPUT=$tempDir/"$barcodeSampleName".bam \
  OUTPUT=$tempDir/"$barcodeSampleName".sorted.bam SORT_ORDER=coordinate


# Add read group and index
java -Xmx4g -jar $picardPath/picard.jar AddOrReplaceReadGroups \
  INPUT=$tempDir/"$barcodeSampleName".sorted.bam \
  OUTPUT=$tempDir/"$barcodeSampleName".sorted.rg.bam \
  RGID=1 RGLB=RGLB RGPL=RGPL RGPU=RGPU RGSM=$barcodeSampleName \
  CREATE_INDEX=True


#mark duplicates
java -Xmx4g -jar $picardPath/picard.jar MarkDuplicates \
      REMOVE_DUPLICATES=true  \
      I=$tempDir/"$barcodeSampleName".sorted.rg.bam \
      O=$tempDir/"$barcodeSampleName".sorted.rg.marked_duplicates.bam \
      M=$tempDir/"$barcodeSampleName"_marked_dup_metrics.txt

#index
java -Xmx4g -jar $picardPath/picard.jar BuildBamIndex \
      INPUT=$tempDir/"$barcodeSampleName".sorted.rg.marked_duplicates.bam


#remove mitochondrial, unmapped and chr Y
samtools view -b $tempDir/"$barcodeSampleName".sorted.rg.marked_duplicates.bam\
 chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX \
 > $tempDir/"$barcodeSampleName".removedChr.bam

# adjust for Tn5 insertion
# forward read: start + 4 , subtract 5 from the partner start
# reverse read: end - 5 , 4 added to its mate start
samtools view -h $tempDir/"$barcodeSampleName".removedChr.bam | \
  awk -v OFS='\t' '{
    if (substr($1,1,1) == "@") { print $0 } else {
      if (and($2, 16) == 0) {$4=$4+4; $9=$9-9};  
      if (and($2, 16) != 0) {$8=$8+4; $9=$9-9};
          print $0 }
  }' | samtools view -b > $tempDir/"$barcodeSampleName".tn5adj.bam

#map quality > 30
samtools view -bq 30 $tempDir/"$barcodeSampleName".tn5adj.bam > \
  $tempDir/"$barcodeSampleName".mapq30.bam

#sort
java -Xmx4g -jar $picardPath/picard.jar SortSam \
  INPUT=$tempDir/"$barcodeSampleName".mapq30.bam \
  OUTPUT=$bamDir/"$barcodeSampleName".final.bam SORT_ORDER=coordinate

#index
samtools index $bamDir/"$barcodeSampleName".final.bam



