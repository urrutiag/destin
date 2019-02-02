destinGrid = function(rse, nClusters,
                      PCrange = 3:25,
                      TSSWeightsList = list(c(1, 2), c(1, 1.5), c(1, 1.25), c(1, 1)),
                      DHSWeightsList = list(c(1, 1), c(1, 2), c(1, 3), c(1, 5)),
                      nCores = NULL, writeOut = F, outDir = NULL){
                      # depthAdjustment = "postPCA" ){
  
  weightGrid = expand.grid(TSSIndex = seq_along(TSSWeightsList), 
                           DHSIndex = seq_along(DHSWeightsList))

  #set up parallel for destin:
  if(!is.null(nCores)){
    cl = makeCluster(nCores)
    clusterEvalQ(cl, library(SummarizedExperiment))
    clusterEvalQ(cl, library(Matrix))
    clusterEvalQ(cl, library(irlba))
    clusterEvalQ(cl, library(data.table))
    clusterExport(cl, list("rse", "TSSWeightsList", "DHSWeightsList", "weightGrid",
                           "getDestin", "PCrange", "getLogLike" ,
                           "nClusters")
                  # , "depthAdjustment")
                  , envir = environment())
    resultsList = parLapply(cl, 1:nrow(weightGrid), function(gridRow) {
      TSSWeights = TSSWeightsList[[weightGrid[gridRow,]$TSSIndex]] 
      DHSWeights = DHSWeightsList[[weightGrid[gridRow,]$DHSIndex]]
      result =  try(getDestin( rse, nClusters = nClusters, PCrange=PCrange, 
                               TSSWeights=TSSWeights, DHSWeights=DHSWeights, )
                              # depthAdjustment = depthAdjustment) 
                    )
      if (class(result) == "try-error") return(NULL)
      return(result)
    } ) 
    stopCluster(cl)
  }
  
  if(is.null(nCores)){
    resultsList = lapply(1:nrow(weightGrid), function(gridRow) {
      TSSWeights = TSSWeightsList[[weightGrid[gridRow,]$TSSIndex]] 
      DHSWeights = DHSWeightsList[[weightGrid[gridRow,]$DHSIndex]]
      result =  try(getDestin( rse, nClusters = nClusters, PCrange=PCrange, 
                               TSSWeights=TSSWeights, DHSWeights=DHSWeights, )
                              # depthAdjustment = depthAdjustment) 
                    )
      if (class(result) == "try-error") return(NULL)
      return(result)
    } ) 
  }
  
  results = rbindlist(resultsList)
  results$sampleName = sampleName
  
  ### Run again with optimal parameters:
  resultsMaxLike = results[which.max(logLike)]
  nPCsOpt = resultsMaxLike$nPCs
  TSSWeightsOpt = c(resultsMaxLike$TSSWeight1, resultsMaxLike$TSSWeight2)
  DHSWeightsOpt = c(resultsMaxLike$DHSWeight1, resultsMaxLike$DHSWeight2)
  resultFinal = try(
    getDestin(rse, PCrange = nPCsOpt, TSSWeight = TSSWeightsOpt, 
              DHSWeight = DHSWeightsOpt, outCluster=T, nClusters = nClusters)
              # depthAdjustment = depthAdjustment)
  )
  
  if(writeOut){
    resultsDir = file.path(outDir, "results")
    write.csv(results, 
              file = file.path(resultsDir, paste0(sampleName, "DestinResultsAllSummary.csv")), 
              row.names = F)

    write.csv(resultFinal$summary, 
              file = file.path(resultsDir, paste0(sampleName, "DestinResultsFinalSummary.csv")), 
              row.names = F)
    
    write.csv(resultFinal$cluster, 
              file = file.path(resultsDir, paste0(sampleName, "DestinResultsFinalCluster.csv")), 
              row.names = F)
  }
  
  return(resultFinal)
} 
