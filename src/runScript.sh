ssh urrutia@longleaf.unc.edu
srun --pty --mem=32000 --time=1000 --cpus-per-task=1 /bin/bash

~/Dropbox/Gene_shared/scATACseq/data
labelsMapP56.txt
p56_cluster.txt

OC 5 - 252
IN2 6 - 329
EX2 7 - 366
# barcodesPre=$(awk '{ print substr($1, 6, 36) }' p56_cluster.txt | head)
#R
cluster = read.table("~/scATACseq/fastqDump/scripts/p56_cluster.txt")
clusterOut = cluster[cluster$V2 %in% c(5,6,7),]
writeLines(text = paste(substr(clusterOut$V1,6,37)),
          con = "~/scATACseq/fastqDump/scripts/barcodes.txt")

#/bin/bash
barcodesPre=$( cat ~/scATACseq/fastqDump/scripts/barcodes.txt )
#barcodesPre=$(awk '{ print }' barcodes.txt | head )
barcodes=( $barcodesPre )



dataDir=~/scATACseq/fastqDump/Preissl
dataDirOG=/pine/scr/u/r/urrutia/fastqDump/Preissl
tempDir=/pine/scr/u/r/urrutia/fastqDump/PreisslTemp

ls $tempDir
rm -r $tempDir
mkdir $tempDir
ls  $dataDir
rm -r $dataDir
mkdir $dataDir

# index1=0; index2=2; myRep=rep1
# sbatch ~/scATACseq/fastqDump/scripts/parseFastq.sh $index1 $index2 $myRep

for index1 in $(seq 0 50 450);  
do    
index2=50;    
if [ $index1 == 450 ]; then index2=35; fi;    
echo $index1;    
echo $index2;    
sbatch ~/scATACseq/fastqDump/scripts/loopThruBarcodes.sh $index1 $index2 rep1;   
done

for index1 in $(seq 485 50 900);  
do    
index2=50;    
if [ $index1 == 900 ]; then index2=38; fi;    
echo $index1;    
echo $index2;    
sbatch ~/scATACseq/fastqDump/scripts/loopThruBarcodes.sh $index1 $index2 rep2;   
done


