
getDestin = function(rse, PCrange=10, TSSWeights=c(1,1), DHSWeights=c(1,1), 
                     nClusters, outCluster = F){
  
  # for normalization
  cellSumPostQC = colData(rse)$cellSumPostQC
  
  #save original accessibility matrix in case region weighting removes regions
  countMatOG = assay(rse)  
  
  ### weight the regions
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
    logLike =  getLogLike(countMatOG, kfit$cluster)
    #if (doLogLikeOG) { logLikeOG =  getLogLike(countMatOG, kfit$cluster) } else {logLikeOG = NA }
    return(list (summary = data.frame(nPCs = myNPC,
                                      #nPeaksNominal = nPeaks,
                                      nPeaksActual = nrow(X),
                                      #pValCutoff = pValCutoff,
                                      purity = purity,
                                      logLike = logLike),
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

### QC ------------------------------------------------------------
doQC = function(rse){
  rse = rse[Matrix::rowSums(assay(rse)) >= 5, ]
  cellSum = Matrix::colSums(assay(rse))
  cutoffs = 2 ^ ( median(log2(cellSum)) + c(-3,3)*mad(log2(cellSum)) )
  rse = rse[ , cellSum > cutoffs[1] & cellSum < cutoffs[2] ]
  cellSumPostQC = Matrix::colSums(assay(rse))
  colData(rse)$cellSumPostQC = cellSumPostQC
  return(rse)
}




