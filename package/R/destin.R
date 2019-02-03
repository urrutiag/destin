
getDestin = function(rse, nClusters, PCrange=10, TSSWeights=c(1,1), DHSWeights=c(1,1),
                     depthAdjustment = "postPCA"){
  
  # for normalization
  cellSumPostQC = colData(rse)$cellSumPostQC
  
  #save original accessibility matrix in case region weighting removes regions
  countMatOG = assay(rse)  
  
  ### weight the regions
  if ( any(TSSWeights != 1 ) | any(DHSWeights != 1 ) ) {
    rowRanges(rse)$TSSMetric[rowRanges(rse)$region == "promoter"] = TSSWeights[1]
    rowRanges(rse)$TSSMetric[rowRanges(rse)$region == "distal element"] = TSSWeights[2]
    rowRanges(rse)$DHSMetric =  dbeta(rowRanges(rse)$DHSsum / 100 + .01, DHSWeights[1], DHSWeights[2])
    rowRanges(rse)$regionWeight =
      (rowRanges(rse)$TSSMetric / mean(rowRanges(rse)$TSSMetric)) *
      (rowRanges(rse)$DHSMetric / mean(rowRanges(rse)$DHSMetric))
    rse = rse[rowRanges(rse)$regionWeight > 0]
  
    X = assay(rse) * rowRanges(rse)$regionWeight  
  } else {
    X = assay(rse)
  }
  
 
  set.seed(10)
  
  ### Adjust for depth pre PCA
  if (depthAdjustment == "prePCA") {
    X_depthPre = t( t(X) / colData(rse)$cellSumPostQC )
    X = X_depthPre
  }

  # PCA   
  pca = irlba(t(X), nv = max(PCrange))
  projection = t(X) %*% pca$v

  ### Adjust for depth post PCA
  if (depthAdjustment == "postPCA") {
    projection = projection / cellSumPostQC
  }

  resultsList = lapply (PCrange, function(myNPC) {
    
    projectionNPCs = projection[, 1:myNPC]
    
    kfit = try(
      kmeans(projectionNPCs, centers = nClusters, nstart = 100)
    )
    if (class(kfit) == "try-error") return (NULL)
    logLike =  getLogLike(countMatOG, kfit$cluster)
    cluster = colData(rse) 
    cluster$cluster = kfit$cluster
    return(list (summary = data.frame(nPCs = myNPC,
                                      TSSWeight1 = TSSWeights[1],
                                      TSSWeight2 = TSSWeights[2],
                                      DHSWeight1 = DHSWeights[1],
                                      DHSWeight2 = DHSWeights[2],
                                      nPeaksActual = nrow(X),
                                     logLike = logLike),
                 cluster = cluster,
                 PCs = projectionNPCs
    ))
  } ) 
  
  return ( resultsList )                          
}


getLogLike = function(countMat, cluster, reassign = FALSE){
  
  cellLikesList = lapply( seq_along( unique(cluster) ), function(myCluster) {
    
    # cluster specific count mat and 
    countMatCluster = countMat[,cluster == myCluster, drop = F]
    
    # empirical probabilities within cluster
    regionTotals = Matrix::rowSums(countMatCluster) 
    if (reassign == T) {
      regionTotals = regionTotals + 1
    }
    regionProbs = regionTotals / sum(regionTotals)
    
    if (reassign == TRUE) {
      countMatEvaluate = countMat
    } else {
      countMatEvaluate = countMatCluster
    }
    
    # term 1: total accessible regions by cell
    cellTotals = colSums(countMatEvaluate)
    term1 = lgamma(cellTotals + 1)
    
    # term 2: sum(log(probabilities)) by cell where probabilities include only accessible regions 
    probsDiag = Diagonal(length(regionProbs), regionProbs)
    # print(regionProbs[1:10])
    scaledCountMat = probsDiag %*% countMatEvaluate
    scaledCountMat@x = log(scaledCountMat@x)
    term2 = colSums(scaledCountMat)
    
    cellLikes = term1 + term2
    return( cellLikes )
    
  })
  
  if (reassign == T){
    logLikeMatrix = do.call(cbind, cellLikesList)
    return( logLikeMatrix )
  }
  
  if (reassign == F){
    clusterLikes = sapply(cellLikesList , sum)
    return( sum ( clusterLikes ) )
  }
  
}



### QC ------------------------------------------------------------

doQC = function(rse, regionSumCutoff = 5, cellSumCutoffSDs = 3){
  rse = rse[Matrix::rowSums(assay(rse)) >= regionSumCutoff, ]
  cellSum = Matrix::colSums(assay(rse))
  cutoffs = 2 ^ ( median(log2(cellSum)) + c(-1,1)*cellSumCutoffSDs*mad(log2(cellSum)) )
  rse = rse[ , cellSum > cutoffs[1] & cellSum < cutoffs[2] ]
  cellSumPostQC = Matrix::colSums(assay(rse))
  colData(rse)$cellSumPostQC = cellSumPostQC
  return(rse)
}



