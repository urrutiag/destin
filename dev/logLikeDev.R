cluster = clusterResult
countMat = assay(rse)[1:10,]

getLogLikeSumFast(countMat, cluster)
getLogLike(countMat, cluster)



getLogLike = function(countMat, cluster, sum = T){
  
  empiricalProbList = lapply(unique(cluster), function(myCluster){
    sums = Matrix::rowSums(countMat[,cluster == myCluster, drop = F]) 
    probs = sums / sum(sums)
    return ( probs )
  })
  names(empiricalProbList) = unique(cluster)
  
  # the col indexing is the bottleneck
  if (sum == T) { 
    logLikes = sapply( seq_along(cluster), function(myCellIndex) {
      dmultFast(x = countMat[,myCellIndex],
                prob = empiricalProbList[[paste(cluster[myCellIndex])]])
    })
    return( sum( logLikes ) )
  } 
  
  # create cell by cluster matrix of likelihoods
  if (sum == F) {
    logLikeList = lapply(1:length(unique(cluster)), function(clusterIndex) {
      logLikes = sapply( seq_along(cluster), function(myCellIndex) {
        dmultFast(x = countMat[,myCellIndex],
                  prob = empiricalProbList[[clusterIndex]])
      })
    })
    logLikeMatrix = do.call(cbind, logLikeList)
    return( logLikeMatrix )
  }
  
}

dmultFast = function(x, prob){ 
  N = sum(x)
  logLike = lgamma(N + 1) + sum(log(prob[x == 1]))
  return( logLike )
}


