destinGrid = function(rse, nClusters,
                      PCrange = 3:25,
                      TSSWeightsList = list(c(1, 2), c(1, 1.5), c(1, 1.25), c(1, 1)),
                      DHSWeightsList = list(c(1, 1), c(1, 2), c(1, 3), c(1, 5)),
                      nCores = NULL, writeOut = F, outDir = NULL, 
                      depthAdjustment = "postPCA" ){
                      
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
                           "getDestin", "PCrange", "getLogLike", "nClusters", 
                           "depthAdjustment")
                  , envir = environment())
    resultsList = parLapply(cl, 1:nrow(weightGrid), function(gridRow) {
      TSSWeights = TSSWeightsList[[weightGrid[gridRow,]$TSSIndex]] 
      DHSWeights = DHSWeightsList[[weightGrid[gridRow,]$DHSIndex]]
      result =  try(getDestin( rse, nClusters = nClusters, PCrange=PCrange, 
                               TSSWeights=TSSWeights, DHSWeights=DHSWeights,
                               depthAdjustment = depthAdjustment) 
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
                               TSSWeights=TSSWeights, DHSWeights=DHSWeights,
                               depthAdjustment = depthAdjustment) 
                    )
      if (class(result) == "try-error") return(NULL)
      return(result)
    } ) 
  }
  
  
  # extract results for highest likelihood
  maxLikelihoodByWeights = sapply(resultsList, function(results) {
    gridSummary = rbindlist( lapply(results , function(result) result$summary ) )
    max(gridSummary$logLike)                     
  } ) 
  resultsMaxWeights = resultsList[[ which.max(maxLikelihoodByWeights) ]]
  
  summarydByPC = rbindlist( lapply(resultsMaxWeights, function(result) result$summary ) )
  resultsFinal =  resultsMaxWeights[[which.max( summarydByPC$logLike ) ]]
  
  
  # bind summaries for all weight / pc combinations
  resultsFinal$summaryAll = rbindlist( lapply( resultsList, function( results ) {
    rbindlist( lapply(results, function(result) result$summary ) )
  } ) ) 
  
  
  if(writeOut){
    resultsDir = file.path(outDir, "results")
    write.csv(resultsFinal$summaryAll, 
              file = file.path(resultsDir, "DestinResultsAllSummary.csv"), 
              row.names = F)

    write.csv(resultFinal$summary, 
              file = file.path(resultsDir, "DestinResultsFinalSummary.csv"), 
              row.names = F)
    
    write.csv(resultFinal$cluster, 
              file = file.path(resultsDir, "DestinResultsFinalCluster.csv"), 
              row.names = F)
  }
  
  return( resultsFinal )
} 
