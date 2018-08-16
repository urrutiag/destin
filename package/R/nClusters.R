
getElbow = function(measureVec, clusterVec){
  elbowMetric  = rep(NA, length(clusterVec))
  for (k in clusterVec) {
    fit1 = lm(measureVec[1:k] ~ clusterVec[1:k])
    fit2 = lm(measureVec[k:length(clusterVec)] ~ clusterVec[k:length(clusterVec)])
    elbowMetric[k] = sum(fit1$residuals^2) + sum(fit2$residuals^2)
  }
  return(elbowMetric)
}

estimateNClusters = function(rse, clusterVec = 1:20){
  
  # catch error, clusterVec must begin with 1
  
  # for normalization
  cellSumPostQC = colData(rse)$cellSumPostQC
  
  #PCA
  PCrange=10
  X = assay(rse)
  pca = irlba(t(X), nv = max(PCrange))
  projection = t(X) %*% pca$v[, 1:PCrange]
  projectionNorm = projection / cellSumPostQC
  projNormMat = as(projectionNorm, "matrix")
  
  #Elbow method:
  logLikeVec = rep(NA, length(clusterVec))
  for (k in clusterVec) { 
    kfit = kmeans(projectionNorm, centers = k)
    logLikeVec[k] = getLogLike(assay(rse), kfit$cluster)
  }
  logLikeElbow = getElbow(logLikeVec, clusterVec)
  
  logLikes = data.frame(nClusters = clusterVec, 
                        logLike = logLikeVec)
  nClustersEstimate =  which.min(logLikeElbow)
  
  #record best fit  need to rework getElbow
  #kfit = kmeans(projectionNorm, centers = nClustersEstimate)
  #logLikeVec[k] = getLogLike(assay(rse), kfit$cluster)
  
  
  
  out = list(logLikes = logLikes,
             nClustersEstimate = nClustersEstimate)
  
  return(out)
}

