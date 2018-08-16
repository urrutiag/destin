#!/usr/bin/env Rscript

#SBATCH --mem=32g
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=24

args = commandArgs(trailingOnly=TRUE)

sampleName = args[1]
outDir = args[2]
cellDataFile = args[3]
model = args[4]

bamDir = file.path(outDir, "bam")
peaksDir = file.path(outDir, "peaks")
bedFile = file.path(peaksDir, 
    paste0(sampleName, "_peaks.blacklist_removed.narrowPeak"))

# source("https://bioconductor.org/biocLite.R")
library(rtracklayer)
library(GenomicAlignments)
library(SummarizedExperiment)
library(Matrix)
library(data.table)
library(ChIPpeakAnno)

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
myBed = import(bedFile, format = "BED",
               extraCols = extraCols_narrowPeak)
elementMetadata(myBed)$name = NULL

  # keep top 500000 regions
  top500quantile = 1 - ( min( 500000, length(myBed) ) / length(myBed) )
  pValCutoff = quantile( myBed$pValue, top200quantile )
  peakBelowCutoff = myBed$pValue > pValCutoff # pvalue on -log scale
  myBed = myBed[peakBelowCutoff,]


bamFiles = dir( bamDir, pattern = "final.bam" )
bamFiles = bamFiles[!grepl("bai", bamFiles)]


#credit to GreenleafLab/chromVAR
left_right_to_grglist <- function(left, right) {
  stopifnot(length(left) == length(right))
  if (length(left) == 0) {
    return(GenomicRangesList())
  }
  x <- c(left, right)[as.vector(matrix(seq_len(2L * length(left)), nrow = 2L, 
                                       byrow = TRUE))]
  p <- PartitioningByEnd(cumsum(rep(2, length(x)/2)))
  out <- relist(x, p)
  return(out)
}


#credit to GreenleafLab/chromVAR
bamToFragments <- function(bamfile, paired) {
    scanned <- scanBam(bamfile, 
                       param = 
                         ScanBamParam(flag = 
                                        scanBamFlag(isMinusStrand = FALSE, 
                                                    isProperPair = TRUE),
                                      what = c("rname", "pos", "isize")))[[1]]
    scanned_left <- GRanges(seqnames = scanned$rname, 
                            IRanges(start = scanned$pos, 
                                    width = 1), strand = "+")
    scanned_right <- GRanges(seqnames = scanned$rname, 
                             IRanges(start = scanned$pos + 
                                       abs(scanned$isize) - 1, width = 1),
                             strand = "-")
    out <- left_right_to_grglist(scanned_left, scanned_right)
  return(out)
}

type = "reads"

countsList = lapply(bamFiles, function( myBamFile ) {
  if (type == "fragments"){
    fragments = bamToFragments(file.path(bamDir, myBamFile), paired = paired)
    counts = as( countOverlaps( myBed, fragments ), "dsparseVector")
  }
  if (type == "reads"){
    bam  = readGAlignments(file.path(bamDir, myBamFile)) 
    counts = as( countOverlaps( myBed, bam ), "dsparseVector")
  }
  list( file = myBamFile, counts = counts)
})

fileNames = sapply(countsList, function(myCounts) {
  myCounts[["file"]]
})

input = lapply(countsList, function(myCounts) {
  myCounts[["counts"]]
})
nRows = length(countsList[[1]]$counts)
countsMat = sparseMatrix( 
            x=unlist(lapply(input,slot,"x")), 
            i=unlist(lapply(input,slot,"i")), 
            p=c(0,cumsum(sapply(input,function(x){length(x@x)}))),
            dims=c(nRows,length(input))
        )
countsMat[countsMat > 1] = 1


cellID = gsub(paste0(".",sampleName,".final.bam"),"", fileNames)
colData = data.table(data.frame(
  cellID = cellID,
  fileName = fileNames
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

rse = SummarizedExperiment(assays=list(counts=countsMat),
                     rowRanges=myBed, colData=colMerged)

metadata(rse) = list(
  dataSet = sampleName
)

rse = rse [ , rse$cell_type != "excluded" ]


### Annotate ------------------------------------------------------------


myFilePath = system.file("annotation/encode", paste0("DHS", model, ".Rdata"), 
                         package="destin")
load(myFilePath)
dhsRanges = DHSobject$genomeRanges
regionIndex = findOverlaps(rowRanges(rse), dhsRanges, select="first")
rowRanges(rse)$DHSsum = DHSobject$sumDHS[regionIndex] 

if ( model == "hg19" ) { 
  data(TSS.human.GRCh37)
  annoFile = TSS.human.GRCh37 
}
if ( model == "mm10" ) {
  data(TSS.mouse.GRCm38)
  annoFile = TSS.mouse.GRCm38  
}
annotatedPeak = annotatePeakInBatch(rowRanges(rse), AnnotationData=annoFile, select = "first" )
rowRanges(rse)$distanceToTSS = annotatedPeak$distancetoFeature
rowRanges(rse)$feature = annotatedPeak$feature
rowRanges(rse)$region = "promoter"
rowRanges(rse)$region[ rowRanges(rse)$distanceToTSS > 1000 | 
                               rowRanges(rse)$distanceToTSS < -2000 ] = "distal element"


save(rse, file = file.path(peaksDir, paste0(sampleName, "500KPeaksAnnotated.Rdata")))



