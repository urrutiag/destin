#!/usr/bin/env Rscript

#SBATCH --mem=32g
#SBATCH --time=24:00:00

args = commandArgs(trailingOnly=TRUE)

sampleName = args[1]
outDir = args[2]
cellDataFile = args[3]
model = args[4]

# sampleName = "BuenrostroMouse"
# outDir = "/proj/yuchaojlab/gene/scATACseq/data/BuenrostroMouse"
# cellDataFile = file.path(outDir, "SraRunTableBuenrostroMouse.txt")
# model = "mm10"

library(BiocParallel)
register(MulticoreParam(1))

# source("https://bioconductor.org/biocLite.R")
# biocLite("remotes")
# biocLite("GreenleafLab/chromVAR")
# biocLite("motifmatchr")
# biocLite("JASPAR2016")

library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2016)
library(SummarizedExperiment)
library(data.table)

bamDir = file.path(outDir, "bam")
peaksDir = file.path(outDir, "peaks")
bedFile = file.path(peaksDir, 
    paste0(sampleName, "_peaks.blacklist_removed.narrowPeak"))


extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
myBed = import(bedFile, format = "BED",
               extraCols = extraCols_narrowPeak)

#keep top 500K regions
nRows = length(myBed)
top500quantile = 1 - ( min( 500000, nRows ) / nRows ) 
pValCutoff = quantile( myBed$pValue, top500quantile ) 
myBed = myBed[ myBed$pValue >= pValCutoff , ]

my.regions <- GRanges(seqnames(myBed),
  IRanges(start(myBed),
  end(myBed)))

bamFiles = dir( bamDir, pattern = "final.bam" )
bamFiles = bamFiles[!grepl("bai", bamFiles)]
bamFiles = bamFiles[!bamFiles %in% 
			c( "SRR3986179.Corces.final.bam",
			"SRR3986211.Corces.final.bam")]

#peaks <- getPeaks(bedFile, 
#  extra_cols=c("signalValue", "pValue",
#                          "qValue", "peak" ))

colData = data.table(data.frame(
  cellID = gsub(paste0(".",sampleName,".final.bam"),"", bamFiles),
  fileName = bamFiles
))

if (cellDataFile != ""){
  cellData = fread(cellDataFile) 
  if ("subSampleID" %in% names(cellData)) {
    colData[, subSampleID := tstrsplit(cellID, ".", fixed = T)[[1]] ]
    cellData[, subSampleID := tstrsplit(bamFile, ".", fixed = T)[[1]] ]
    colMerged = merge(colData, cellData, by = "subSampleID")
  } else {
    colMerged = merge(colData, cellData, by.x = "cellID", by.y = "Run")
  }
} else {
  colMerged = colData
}

#sparse Matrix of class "dgCMatrix"
setwd(bamDir)
rse <- getCounts(bamFiles, my.regions, 
                              paired =  TRUE, 
                              by_rg = F, 
                              format = "bam", 
                              colData = colMerged)
assay(rse)[assay(rse) > 1] = 1

if ( model == "hg19" ) {
  myGenome = BSgenome.Hsapiens.UCSC.hg19
  mySpecies = "Homo sapiens"
}
if ( model == "mm10" ) {
  myGenome = BSgenome.Mmusculus.UCSC.mm10
  mySpecies = "Mus musculus"
}

rse <- addGCBias(rse, genome = myGenome)
rse_filtered <- filterSamples(rse, min_depth = 1500,
                                  min_in_peaks = 0.15, shiny = F)
rse_filtered = sort(rse_filtered)
rse_filtered <- filterPeaks(rse_filtered)
motifs <- getJasparMotifs(species = mySpecies)
motif_ix <- matchMotifs(motifs, rse_filtered,
                         genome = myGenome)

# computing deviations
dev <- computeDeviations(object = rse_filtered, 
                                 annotations = motif_ix)


classes = colMerged$cell_type[
  match(colData(dev)$fileName, colMerged$fileName)]
nClusters = length(unique(classes))

#hclust
sample_dist <- getSampleDistance(dev) 
clustResult = hclust(sample_dist)
memb <- cutree(clustResult, k = nClusters)
hTable = table(memb, classes)
purityHclust = sum(apply(hTable,1,max))/sum(hTable)

#kmeans
kfit = kmeans(t(assay(dev)), centers = nClusters, nstart=100)
kTable = table(kfit$cluster, classes)
purityKmeans = sum(apply(kTable,1,max))/sum(kTable)

out = data.frame(fileName = colData(dev)$fileName,
                 cellID = colData(dev)$cellID,
                 cell_type = colData(dev)$cell_type,
                 clusterHclust = memb,
                 clusterKmeans = kfit$cluster,
                 purityHclust = purityHclust,
                 purityKmeans = purityKmeans
)

write.csv(out, file = file.path(outDir, "peaks", 
  paste0(sampleName, "ChromVarAllPeaks.csv")))

