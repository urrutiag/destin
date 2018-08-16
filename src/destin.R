#!/usr/bin/env Rscript

#SBATCH --mem=32g
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=24

args = commandArgs(trailingOnly=TRUE)

sampleName = args[1]
outDir = args[2]
cellDataFile = args[3]
# model = args[4]

# sampleName = "BuenrostroHuman"
# outDir = "/proj/yuchaojlab/gene/scATACseq/data/BuenrostroHuman"
# cellDataFile = file.path(outDir, "SraRunTableBuenrostroHuman.txt")
# model = "hg19"
# gridRow = 1

peaksDir = file.path(outDir, "peaks")

library(data.table)
library(ggplot2)
library(irlba)
library(SummarizedExperiment)
library(rtracklayer)
library(ChIPpeakAnno)
library(parallel)

getDestin = function(rse, PCrange=10, TSSWeights=c(1,1), DHSWeights=c(1,1), 
                     nPeaks=1e5, pValCutoff=NULL, nClusters=NULL, outCluster = F, doLogLikeOG=F){

  countMatOG = assay(rse)  
  if (is.null(nClusters)) nClusters = length(unique(rse$cell_type))

  ### keep top nPeaks regions
  if (is.null( pValCutoff)){
    top100quantile = 1 - ( min( nPeaks, nrow(rse) ) / nrow(rse) )
    pValCutoff = quantile( rowRanges(rse)$pValue, top100quantile )
  } 
  peakBelowCutoff = rowRanges(rse)$pValue >= pValCutoff # pvalue on -log scale
  if (!is.null( pValCutoff))  nPeaks = sum(peakBelowCutoff)
  rse = rse[peakBelowCutoff,]
  
  ### QC ------------------------------------------------------------
  rse = rse[Matrix::rowSums(assay(rse)) >= 5, ]
  cellSum = Matrix::colSums(assay(rse))
  cutoffs = 2 ^ ( median(log2(cellSum)) + c(-3,3)*mad(log2(cellSum)) )
  rse = rse[ , cellSum > cutoffs[1] & cellSum < cutoffs[2] ]
  cellSumPostQC = Matrix::colSums(assay(rse))

  ### weight the loci
  rowRanges(rse)$TSSMetric[rowRanges(rse)$region == "promoter"] = TSSWeights[1] 
  rowRanges(rse)$TSSMetric[rowRanges(rse)$region == "distal element"] = TSSWeights[2]
  rowRanges(rse)$DHSMetric =  dbeta( rowRanges(rse)$DHSsum/100 + .01, DHSWeights[1], DHSWeights[2]) 
  rowRanges(rse)$regionWeight = 
    (rowRanges(rse)$TSSMetric / mean(rowRanges(rse)$TSSMetric)) * 
    (rowRanges(rse)$DHSMetric / mean(rowRanges(rse)$DHSMetric)) 
  rse = rse[rowRanges(rse)$regionWeight > 0]
  
  
  # DR and cluster ---------------------------------------------
  set.seed(10)

  X = assay(rse) * rowRanges(rse)$regionWeight  
  pca = try(
    irlba(t(X), nv = max(PCrange))
  )
  if (class(pca) == "try-error") return (NULL)
  
  results = lapply (PCrange, function(myNPC) {
    projection = t(X) %*% pca$v[, 1:myNPC]
    projectionNorm = projection / cellSumPostQC
    kfit = try(
      kmeans(projectionNorm, centers = nClusters, nstart = 100)
    )
    if (class(kfit) == "try-error") return (NULL)
    purity = sum( apply(table(kfit$cluster, rse$cell_type), 1, max) ) / ncol(rse)
    logLike =  getLogLike(assay(rse), kfit$cluster)
    if (doLogLikeOG) { logLikeOG =  getLogLike(countMatOG, kfit$cluster) } else {logLikeOG = NA }
    return(list (summary = data.frame(nPCs = myNPC,
	  	      nPeaksNominal = nPeaks,
                      nPeaksActual = nrow(X),
                      pValCutoff = pValCutoff,
                      purity = purity,
                      logLike = logLike,
                      logLikeOG = logLikeOG),
          cluster = data.frame(cellID = rse$cellID,
                               cluster = kfit$cluster)
          ))
  } ) 
  
  resultsSummary = rbindlist(lapply(results, function(x) x$summary))

  summary = data.frame(
    TSSWeight1 = TSSWeights[1],
    TSSWeight2 = TSSWeights[2],
    DHSWeight1 = DHSWeights[1],
    DHSWeight2 = DHSWeights[2],
    resultsSummary)
  
  if (outCluster == F) return(summary)

  if (outCluster == T) return(list(summary = summary, cluster = results[[1]]$cluster))
}


getLogLike = function(countMat, cluster){
  
  empiricalProbList = lapply(unique(cluster), function(myCellType){
    sums = Matrix::rowSums(countMat[,cluster == myCellType, drop = F]) 
    probs = sums / sum(sums)
    return ( probs )
  })
  names(empiricalProbList) = unique(cluster)
  
  logLikes = sapply( seq_along(cluster), function(myCellIndex) {
    dmultFast(x = countMat[,myCellIndex],
              prob = empiricalProbList[[paste(cluster[myCellIndex])]])
  })
  
  return(sum(logLikes))
}

dmultFast = function(x, prob){ 
  N = sum(x)
  logLike = lgamma(N + 1) + sum(log(prob[x == 1]))
  return( logLike )
}




#begin---------------------------------------------------------------------

#load data:
fpath <- file.path(peaksDir, paste0(sampleName, "500KPeaksAnnotated.Rdata"))
load(fpath)


runDestinOG = T
# Begin Destin OG
if (runDestinOG) {
#Hyperparameters:
PCrange = 3:25
TSSWeightsList = list(c(1,2), c(1,1.5), c(1,1.25), c(1,1))
DHSWeightsList = list(c(1,1), c(1,2), c(1,3), c(1,5))
weightGrid = expand.grid(TSSIndex = seq_along(TSSWeightsList), 
            DHSIndex = seq_along(DHSWeightsList))
#nPeaksVec = c(5e4, 10e4, 15e4, 20e4)[2]

#set up parallel for destin:
nCores = 20
cl = makeCluster(nCores)
clusterEvalQ(cl, library(SummarizedExperiment))
clusterEvalQ(cl, library(Matrix))
clusterEvalQ(cl, library(irlba))
clusterEvalQ(cl, library(data.table))
clusterExport(cl, list("rse", "TSSWeightsList", "DHSWeightsList", "weightGrid",
                      "getDestin", "PCrange", "getLogLike","dmultFast" ,
                      "peaksDir", "sampleName"))

### This is the current version of Destin
resultsList = parLapply(cl, 1:nrow(weightGrid), function(gridRow) {
   TSSWeights = TSSWeightsList[[weightGrid[gridRow,]$TSSIndex]] 
   DHSWeights = DHSWeightsList[[weightGrid[gridRow,]$DHSIndex]]
#   result = try( rbindlist(  
#     lapply(nPeaksVec, function(nPeaks) {
   result =  try(getDestin(rse, PCrange=PCrange, TSSWeights=TSSWeights, DHSWeights=DHSWeights, 
                 pValCutoff = .01))# , 
#              nPeaks = nPeaks,  doLogLikeOG=T)
#   } ) ) )
  if (class(result) == "try-error") return(NULL)
  result$sampleName = sampleName
  return(result[])
} ) 
results = rbindlist(resultsList)

stopCluster(cl)

write.csv(results, 
         file = file.path(peaksDir, paste0(sampleName, "DestinResultsAllSummary.csv")), 
         row.names = F)

### Run again with optimal parameters:
resultsMaxLike = results[which.max(logLike)]
nPCsOpt = resultsMaxLike$nPCs
#nPeaksOpt = resultsMaxLike$nPeaksNominal
TSSWeightsOpt = c(resultsMaxLike$TSSWeight1, resultsMaxLike$TSSWeight2)
DHSWeightsOpt = c(resultsMaxLike$DHSWeight1, resultsMaxLike$DHSWeight2)

resultFinal = try(
     getDestin(rse, PCrange = nPCsOpt, TSSWeight = TSSWeightsOpt, 
               DHSWeight = DHSWeightsOpt, outCluster=T, pValCutoff = 2)
   )

write.csv(resultFinal$summary, 
 file.path(peaksDir, paste0(sampleName, "DestinResultsFinalSummary.csv")), 
 row.names = F)

write.csv(resultFinal$cluster, 
  file.path(peaksDir, paste0(sampleName, "DestinResultsFinalCluster.csv")), 
  row.names = F)

} #end Destin OG




###Diagnostics-----------------------------------------

nPeaksVec = seq(2e4,nrow(rse),2e4)
pValCutoffVec = 2:40
TSSWeightsList = list(c(0,1), c(1,4), c(1,2), c(1,1.5), c(1,1.25), c(1,1),
                      c(1.25,1), c(1.5,1), c(2,1), c(4,1), c(1,0))
DHSWeightsList = list(c(5,1), c(4,1), c(3,1), c(2,1), c(1,1), c(1,2), c(1,3),
                      c(1,4), c(1,5))
PCrange = 3:50

nCores = 24
cl = makeCluster(nCores)
clusterEvalQ(cl, library(SummarizedExperiment))
clusterEvalQ(cl, library(Matrix))
clusterEvalQ(cl, library(irlba))
clusterEvalQ(cl, library(data.table))
clusterExport(cl, list("rse", "TSSWeightsList", "DHSWeightsList", 
                      "getDestin", "PCrange", "getLogLike", "dmultFast",
                      "peaksDir", "sampleName", "nPeaksVec", "pValCutoffVec"))
 
pValCutoffResults  =  rbindlist(
  parLapply(cl, pValCutoffVec, function(pValCutoff) {
        getDestin(rse, pValCutoff = pValCutoff, doLogLikeOG=T)
}))
pValCutoffResults$type = "pValCutoff"
nPeaksResults  =  rbindlist(
  parLapply(cl, nPeaksVec, function(myPeak) {
        getDestin(rse, nPeaks = myPeak, doLogLikeOG=T)
}))
nPeaksResults$type = "nPeaks"
TSSWeightsResults  =  rbindlist(
  parLapply(cl, TSSWeightsList, function(TSSWeights) {
        getDestin(rse, TSSWeights = TSSWeights, doLogLikeOG=T)
}))
TSSWeightsResults$type = "TSSWeights"
DHSWeightsResults  =  rbindlist(
  parLapply(cl, DHSWeightsList, function(DHSWeights) {
        getDestin(rse, DHSWeights = DHSWeights, doLogLikeOG=T)
}))
DHSWeightsResults$type = "DHSWeights"
nPCResults =  rbindlist(
  parLapply(cl, PCrange, function(nPCs) {
        getDestin(rse, PCrange = nPCs, doLogLikeOG=T)
}))
nPCResults$type = "nPCs"

stopCluster(cl)

diagnosticsResults = rbind(
  pValCutoffResults,
  nPeaksResults,
  TSSWeightsResults,
  DHSWeightsResults,
  nPCResults
)

write.csv(diagnosticsResults,
 file.path(peaksDir, paste0(sampleName, "DestinDiagnostics.csv")),
 row.names = F)

optList = list(
  nPeaks = nPeaksResults[which.max(logLikeOG),]$nPeaksNominal,
  TSSWeights = c(TSSWeightsResults[which.max(logLikeOG),]$TSSWeight1,
                 TSSWeightsResults[which.max(logLikeOG),]$TSSWeight2),
  DHSWeights = c(DHSWeightsResults[which.max(logLikeOG),]$DHSWeight1,
                 DHSWeightsResults[which.max(logLikeOG),]$DHSWeight2),
  nPCs = nPCResults[which.max(logLikeOG),]$nPCs
)

resultFinalLinear = try(
     getDestin(rse, PCrange = optList$nPCs, TSSWeight = optList$TSSWeights,
               DHSWeight = optList$DHSWeights, nPeaks = optList$nPeaks)
   )

write.csv(resultFinalLinear,
 file.path(peaksDir, paste0(sampleName, "DestinResultsLinearSummary.csv")),
 row.names = F)

