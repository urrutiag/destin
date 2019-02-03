subset 

# nCells = 100
# nRegions = 1000
# 
# # subset cells
# cellSamples = sample(1:ncol(rse), nCells)
# rse = rse[,cellSamples]
# 
# # subset regions
# regionSamples = sample(1:nrow(rse), nRegions)
# rse = rse[regionSamples,]



```{r}
tempDir = '~/Documents/temp'

destinGridBig = function(rse, sampleName, tempDir,
                         PCrange = 3:25,
                         TSSParamList = list(c(1,2), c(1,1.5), c(1,1.25), c(1,1)),
                         DHSParamList = list(c(1,1), c(1,2), c(1,3), c(1,5)),
                         nClusters, nCores=NULL, writeOut=F, outDir=NULL){
  # write count mat
  dir.create(tempDir, showWarnings = F)
  Matrix::writeMM(assay(rse), file = file.path(tempDir, 'countMat.mtx'))
  
  cellSumPostQC = colData(rse)$cellSumPostQC
  
  #calculate and write weights
  weightGrid = expand.grid(TSSIndex = seq_along(TSSWeightsList), 
                           DHSIndex = seq_along(DHSWeightsList))
  weightMatrix = sapply(1:nrow(weightGrid), function(gridRow) {
    TSSparams = TSSParamList[[weightGrid[gridRow,]$TSSIndex]] 
    DHSparams = DHSParamList[[weightGrid[gridRow,]$DHSIndex]]
    TSSweights = rep(0, nrow(rse))
    TSSweights[rowRanges(rse)$region == "promoter"] = TSSparams[1]
    TSSweights[rowRanges(rse)$region == "distal element"] = TSSparams[2]
    DHSweights =  dbeta( rowRanges(rse)$DHSsum/100 + .01, DHSparams[1], DHSparams[2])
    # rse = rse[rowRanges(rse)$regionWeight > 0]
    regionWeights = ( TSSweights / mean(TSSweights) ) * ( DHSweights / mean(DHSweights) ) 
    return(regionWeights)
  } )
  write.csv(weightMatrix, file = file.path(tempDir, 'weightMatrix.mtx'), row.names = F)
}