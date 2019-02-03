# TODO: add purity

yourPathToDestinRepo = "~/Documents/gitRepos/destin"
install.packages(file.path(yourPathToDestinRepo,"package"), repos = NULL, type = "source")

library(destin, quietly = T)

library(proxy)

library(scABC)

library(BiocParallel)
register(MulticoreParam(7))

library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2016)

functionPath = file.path(yourPathToDestinRepo, "archive/timeFunctions.R")
source(functionPath)

timeitSingleDataset = function(dataIndex, nCells = NULL) {

  sampleNames = c('BuenrostroHuman', 'BuenrostroMouse', 'Corces', 'p56')
  pipeResultsDirs = file.path('~/Documents/bamData', sampleNames)
  models = c('hg19', 'mm10', 'hg19', 'mm10')
  nClustersList = c(6, 2, 4, 8)
  
  pipeResultsDir = pipeResultsDirs[dataIndex]
  sampleName = sampleNames[dataIndex]
  model = models[dataIndex]
  nClusters = nClustersList[dataIndex]
  
  # cellData = NULL
  # if (sampleName == 'BuenrostroMouse') {
  #   cellDataFile = '~/Dropbox/Documents/uncDesktop/sraRuns/SraRunTableBuenrostroMouse.txt'
  #   cellData = fread(cellDataFile)
  #   cellData = cellData[, c('Run', 'source_name')]
  # }
  
  algorithms = c("Destin", "Scasat", "scABC", "chromVar")[1]
  resultList = lapply(algorithms, function(algorithm) {
    # try( 
      runAlgorithm(pipeResultsDir,
                      sampleName,
                      model,
                      nClusters,
                      nCells,
                      algorithm) 
    # )
  })
  results = rbindlist(resultList)
  
  return(results)
}

nCells = NULL
resultsList = list()
for (dataIndex in 1:4) {
  resultsList[[dataIndex]] = timeitSingleDataset(dataIndex, nCells = nCells)
}


# 10x datasets

library(destin)
functionPath = "~/Dropbox/Gene_shared/scATACseq/scripts/timeFunctions.R"
source(functionPath)

sampleNames = c('atac_v1_pbmc_5k', 'atac_v1_pbmc_10k')
nCells = NULL
nRegions = NULL
for (sampleName in sampleNames) {
  cluster10x(sampleName, nCells = nCells, nRegions = nRegions)
}


#-------------- Cusanovich 2018

library(data.table)
library(SummarizedExperiment)

dataDir = "~/Documents/Cusanovich2018/Cerebellum"
fileNames = dir(dataDir, full.names = T)
windowmatrixName = fileNames[grepl('5kbwindowmatrix', fileNames)]

# indextableName = fileNames[grepl('indextable', fileNames)]
# indexTable = fread(indextableName, header = F)
# setnames(indexTable, c('barcode', 'sample'))

chunksize = 100000 # 540470 regions total
i = 1
nskip = (i-1)*chunksize

windowmatrix = fread(windowmatrixName, nrow = chunksize, skip = nskip)
format( object.size(windowmatrix), units="Mb" )

countmat = as.matrix(windowmatrix[,6:ncol(windowmatrix)])
countMat = as(countmat, 'sparseMatrix')
format( object.size(countMat), units="Mb")



# 44 Mb for 100K regions, 3K cells
# 9 Gb for 600K regions, 100K cells

can possibly read into python 

cellID = colnames(windowmatrix)[6:ncol(windowmatrix)]





windowmatrixSub = windowmatrix[1:1000,1:1000]
# indexTableSub = indexTable[1:(10 - 5)]

rowRanges = GRanges(windowmatrixSub[, c('chr', 'start',   'end')])

windowmatrixSub[,c('V1', 'chr', 'start', 'end', 'annot') := NULL]
colData = data.frame(cellID = names(windowmatrixSub))
countsMat = as(windowmatrixSub, 'matrix')
counts <- as(countsMat , "sparseMatrix") 

rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
                            rowRanges=rowRanges, colData=colData)



