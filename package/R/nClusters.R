
getElbow = function(measureVec, clusterVec){
  elbowMetric  = rep(NA, length(clusterVec))
  for (k in clusterVec) {
    fit1 = lm(measureVec[1:k] ~ clusterVec[1:k])
    fit2 = lm(measureVec[k:length(clusterVec)] ~ clusterVec[k:length(clusterVec)])
    elbowMetric[k] = sum(fit1$residuals^2) + sum(fit2$residuals^2)
  }
  return(elbowMetric)
}

#TODO check error on range starts with 2
estimateNClusters = function(rse, nClustersRange = 2:20, allMethods = F){
  
  clusterVec = nClustersRange
  
  # catch error, clusterVec must begin with 1
  # if (clusterVec[1] != 1){
  #   stop("Please ensure that clusterVec begins with 1")
  # }
  
  # for normalization
  cellSumPostQC = colData(rse)$cellSumPostQC

  #PCA
  PCrange=10
  
  X = assay(rse)
  
  pca = irlba(t(X), nv = max(PCrange))
  
  projection = t(X) %*% pca$v[, 1:PCrange]
  projectionNorm = projection / cellSumPostQC
  projNormMat = as(projectionNorm, "matrix")
  
  metricsList = list()
  nClustersList = list()
  
  # model-based likelihood elbow
  # WCSSE elbow
  print('computing elbow methods')
  logLikeVec = rep(NA, length(clusterVec))
  wcsseVec = rep(NA, length(clusterVec))
  for (k in clusterVec) { 
    kfit = kmeans(projectionNorm, centers = k)
    logLikeVec[k] = getLogLike(assay(rse), kfit$cluster)
    wcsseVec[k] = kfit$tot.withinss
  }
  metricsList$logLikeElbow = getElbow(logLikeVec, clusterVec)
  metricsList$wcsseElbow = getElbow(wcsseVec, clusterVec)
  nClustersList$logLikeElbow = which.min(metricsList$logLikeElbow)
  nClustersList$wcsseElbow = which.min(metricsList$wcsseElbow)
  
  if (allMethods == T){
    
    #  silhouette
    print('computing silhouette')
    myDist = dist(projNormMat)
    silhouetteStatsNot1 = sapply (2:max(clusterVec), function(myK)
      { 
      clusterKmeans = kmeans(projNormMat, myK, nstart = 100)
      silhouetteOut = cluster::silhouette(clusterKmeans$cluster, myDist)
      return ( summary(silhouetteOut)$avg.width )
    })
    metricsList$silhouetteStats = c(NA, silhouetteStatsNot1)
    nClustersList$silhouette = which.max(metricsList$silhouetteStats)
    
    #distortion
    if ( "ClusterR" %in% rownames(installed.packages()) ) {
      print('computing distortion')
      distortionStats = try(
        ClusterR::Optimal_Clusters_KMeans(
          data = projNormMat, 
          max_clusters = max(clusterVec), 
          plot_clusters = F, 
          criterion = "distortion_fK"),
          silent = T)
      if (class(class(metricsList$distortionStats)) != "try-error")
      {
        metricsList$distortionStats = as.vector(distortionStats)
        nClustersList$distortion = which.min(metricsList$distortionStats)
      } 
    } else {
      print("ClusterR package not installed: distortion statistic unavailable")
    }
    

      #gap
    
      print('computing gap')
      gap_statKmeans = try(
        cluster::clusGap(projNormMat, FUN = kmeans,
                                        K.max = max(clusterVec), B = 500)
        , silent = T)
      if (class(class(metricsList$distortionStats)) != "try-error")
      {
        metricsList$gap_statKmeans = gap_statKmeans$Tab[,3]
        nClustersList$GapStat = cluster::maxSE(gap_statKmeans$Tab[,3], 
                                               gap_statKmeans$Tab[,4], 
                                               method = "firstSEmax")
      }
    
  }
  
  out = list(nClustersList = nClustersList,
             metricsList = metricsList)
  
  return(out)
}


plotNClusters = function(clusterEst){
  nClustersMax = length(clusterEst$metricsList[[1]])
  pList = list()
  for (j in seq_len(length(clusterEst$metricsList))){
    select = rep(FALSE, nClustersMax)
    select[clusterEst$nClustersList[[j]]] = TRUE
    myData = data.frame(x = 1:nClustersMax, 
                        y = clusterEst$metricsList[[j]], 
                        select = select)
    pList[[j]] = qplot(data = myData, x = x, y = y, color = select,
                       ylab = names(clusterEst$nClustersList)[j], 
                       xlab = "N Clusters") +
      guides(color = FALSE)
  }
  grid.arrange(grobs = pList, ncol = 2)
}


plotCluster = function(clusterResults, type = "PCA", clusterLabels = NULL){

  clusterAssignment = clusterResults$cluster$cluster
  PCs = clusterResults$PCs
  
  if (type == "t-SNE") {
    tsne = Rtsne( as.matrix(PCs) )
    componentX = tsne$Y[,1]
    componentY = tsne$Y[,2]
  }
  if (type == "PCA") {
    componentX = PCs[,1]
    componentY = PCs[,2]
  }
  
  p = qplot(x = componentX,
            y = componentY,
            col = factor(clusterAssignment),
            xlab = paste(type, "component 1"),
            ylab = paste(type, "component 2"),
            main = paste0( "2D ", type, " Visualization")) +  
    guides(color=FALSE)
  
  #annotation optional
  if ( !is.null( clusterLabels ) ){
    clusterCenters = data.frame(
      label = clusterLabels,
      x = tapply(componentX, clusterAssignment, mean),
      y = tapply(componentY, clusterAssignment, mean)
    )
    p = p + with(clusterCenters, annotate("text", x=x, y=y, label = label, size = 6))
  }
  
  print(p)
}
