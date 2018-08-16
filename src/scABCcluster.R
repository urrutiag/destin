#!/usr/bin/env Rscript

#SBATCH --mem=32g
#SBATCH --time=24:00:00

library(scABC)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

sampleName = args[1]
outDir = args[2]
cellDataFile = args[3]

#sampleName = "BuenrostroHuman"
#outDir = "/proj/yuchaojlab/gene/scATACseq/data/BuenrostroHuman"
#cellDataFile = "/proj/yuchaojlab/gene/scATACseq/data/BuenrostroHuman/SraRunTableBuenrostroHuman.txt"

bamDir = file.path(outDir, "bam")
peaksDir = file.path(outDir, "peaks")

if (length(args) == 2){
  cellDataFile = ""
}

bamFiles = dir( bamDir, pattern = "final.bam" )
bamFiles = bamFiles[!grepl("bai", bamFiles)]
bedFile = file.path(peaksDir, 
    paste0(sampleName, "_peaks.blacklist_removed.narrowPeak"))

cellData = fread(cellDataFile)
keepIndex = cellData$cell_type != "excluded"
nClusters = length(unique(cellData$cell_type[keepIndex]))
if ( grepl ( "p56", sampleName ) ) {
  bamFilesKeep = paste0(cellData$Run[keepIndex], ".", sampleName, ".final.bam")
  bamFiles = bamFiles[bamFiles %in% bamFilesKeep]
}


setwd(bamDir)
# s1 = scABC(bamFiles, bedFile, nClusters = 2:15)
# clusterEstimate = length(unique(s1$cluster_assignments))

setwd(bamDir)
s2 = scABC(bamFiles, bedFile, nClusters = nClusters)

cellID = gsub(paste0(".",sampleName,".final.bam"),"", names(s2$cluster_assignments))

out = data.frame(cluster = s2$cluster_assignments,
                 cellID = cellID)

if (cellDataFile != "") {

  if ("subSampleID" %in% names(cellData)){
    cellType = cellData$cell_type[match( names(s2$cluster_assignments), cellData$bamFile )]
  } else {
    cellType = cellData$cell_type[match( cellID, cellData$Run )]
  }
 
  clusterTable = table(s2$cluster_assignments, cellType)
  purity = sum(apply(clusterTable, 1, max))/sum(clusterTable)

  out = data.frame(out,
                   cellType = cellType,
                   purity = purity,
                   nDropped = length(bamFiles) - length(s2$cluster_assignments)
                  #  clusterEstimate = clusterEstimate
  )
}

write.csv(out, file = file.path(peaksDir, paste0(sampleName, "scABCresultsAllPeaks.csv")),
          row.names = F, col.names = F)


