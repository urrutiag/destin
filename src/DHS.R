
dataDir = "~/Documents/encodeDump/bigDump/encodeDump"

# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")
# biocLite("BSgenome.Mmusculus.UCSC.mm10")
library(rtracklayer)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
library(parallel)

#create template bin 500bp
model = "mm10"
if (model == "hg19"){
  genome <- BSgenome.Hsapiens.UCSC.hg19
  nChrom = 23
  dataDir = "~/Documents/encodeDump/bigDump/encodeDump"
}
if (model == "mm10"){
  genome <- BSgenome.Mmusculus.UCSC.mm10
  nChrom = 20
  dataDir = "~/Documents/encodeDump/bigDump/dumpMM10"
} 

gr <- GRanges(
  seqnames=seqnames(genome)[1:nChrom],
  ranges=IRanges(1, end = seqlengths(genome)[1:nChrom]),
  strand="*",
  seqlengths=seqlengths(genome)[1:nChrom])
genomeTiled = unlist(tile(gr, width = 500))


fileNames = dir(dataDir, full.names=T)
myBeds = fileNames[ !grepl("tsv", fileNames)]
File.accession = substr(myBeds, 54, 64)
fileNameMap = data.table(data.frame(
  filePath = myBeds, 
  File.accession = File.accession
) )
myTSV = fileNames[ grepl("tsv", fileNames)]
tsvData = data.table( read.table(myTSV, fill = T, header = T, sep = "\t") )
tsvData = merge( tsvData, fileNameMap, by = "File.accession" )


cellTypes = unique(tsvData$Biosample.term.name)
extraCols_broadPeak = c(signalValue = "numeric", pValue = "numeric",
                        qValue = "numeric")

nCores = detectCores() - 1
cl = makeCluster(nCores,  type="PSOCK")
clusterEvalQ(cl, {
  library(rtracklayer)
  library(data.table)
  library(BSgenome.Hsapiens.UCSC.hg19)
})
clusterExport(cl, list("cellTypes", "tsvData", 
                       "genomeTiled", "extraCols_broadPeak"))

coverageList =  parLapply(cl, cellTypes, function( myCell ){ 
  cellFilePaths = tsvData[Biosample.term.name == myCell]$filePath

  overlapList = lapply(cellFilePaths, function(myCellFilePath) {
    dhsRanges = rtracklayer::import(paste(myCellFilePath), format = "BED",
                                    extraCols = extraCols_broadPeak)
   myOverlap = overlapsAny(genomeTiled, dhsRanges)
   return( myOverlap )
  })

  out = list( meanCoverage = apply( do.call( cbind, overlapList ) ,1, mean ) )
  rm(overlapList)
  names(out) = myCell
  return(out) 
} )

stopCluster(cl)

library(Matrix)
myData = elementMetadata(genomeTiled)= as.data.frame(coverageList)
myDataSparse = Matrix(as.matrix(myData), sparse = T)

# prac  = elementMetadata(genomeTiled)[1:1000,1]
# prac.sparse = Matrix(prac, sparse = T)

class(myData)
print(object.size(myData), units="auto")
class(myDataSparse)
print(object.size(myDataSparse), units="auto")

genomeRangesNull = genomeTiled
elementMetadata(genomeRangesNull) = NULL

sumDHS = apply(myData, 1, sum)
# sumSparse = sumDHSMatrix(sumDHS, sparse = T) not helpful

DHSobject = list(
  metadata = list(
    type = "tissue", 
    binSize = 500,
    source = "ENCODE",
    fileType = "broadPeak Hotspot",
    model = "mm10",
    assay = "DNase-seq"
  ),
  tsv = tsvData,
  genomeRanges = genomeRangesNull,
  DHSmeanCoverage = myDataSparse,
  sumDHS = sumDHS
)

DHSranges =  DHSobject$genomeRanges
elementMetadata(DHSranges)$DHSfrequency = DHSobject$sumDHS
metadata(DHSranges) = DHSobject$metadata
metadata(DHSranges)$tsv = DHSobject$tsv

save(DHSranges, file = "~/Documents/gitRepos/destin/package/inst/annotation/encode/DHShg19.Rdata")
