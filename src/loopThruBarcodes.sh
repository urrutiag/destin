#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH --mem=32000

index1=$1
index2=$2
myRep=$3

barcodes=( $( cat ~/scATACseq/fastqDump/scripts/barcodes.txt ) )
barcodesSub=( ${barcodes[@]:$index1:$index2} )

for barcode in ${barcodesSub[@]}
 do
 echo $barcode
 echo $myRep
   ~/scATACseq/fastqDump/scripts/codePreissl.sh $barcode $myRep
 done

