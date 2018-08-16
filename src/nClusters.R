#!/usr/bin/env Rscript

#SBATCH --mem=32g
#SBATCH --time=24:00:00

args = commandArgs(trailingOnly=TRUE)

sampleName = args[1]
outDir = args[2]

peaksDir = file.path(outDir, "peaks")

library(data.table)
library(irlba)
library(SummarizedExperiment)
library(ClusterR)
library(cluster)

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


getElbow = function(measureVec, clusterVec){
  elbowMetric  = rep(NA, length(clusterVec))
  for (k in clusterVec) {
    fit1 = lm(measureVec[1:k] ~ clusterVec[1:k])
    fit2 = lm(measureVec[k:length(clusterVec)] ~ clusterVec[k:length(clusterVec)])
    elbowMetric[k] = sum(fit1$residuals^2) + sum(fit2$residuals^2)
  }
  return(elbowMetric)
}


fpath <- file.path(peaksDir, paste0(sampleName, "500KPeaksAnnotated.Rdata"))
load(fpath)

PCrange=30

### QC ------------------------------------------------------------
rse = rse[Matrix::rowSums(assay(rse)) >= 5, ]
cellSum = Matrix::colSums(assay(rse))
cutoffs = 2 ^ ( median(log2(cellSum)) + c(-3,3)*mad(log2(cellSum)) )
rse = rse[ , cellSum > cutoffs[1] & cellSum < cutoffs[2] ]
cellSumPostQC = Matrix::colSums(assay(rse))

X = assay(rse)
pca = irlba(t(X), nv = max(PCrange))
projection = t(X) %*% pca$v[, 1:PCrange]
projectionNorm = projection / cellSumPostQC
projNormMat = as(projectionNorm, "matrix")

nClustDistortion = ClusterR::Optimal_Clusters_KMeans(
  projNormMat, max_clusters = 20, plot_clusters = F,
  criterion = 'distortion_fK')
nClustSilhouette = ClusterR::Optimal_Clusters_KMeans(
  projNormMat, max_clusters = 20, plot_clusters = F,
  criterion = 'silhouette')
gap_statKmeans = cluster::clusGap(projNormMat, FUN = kmeans,
                                    K.max = 20, B = 500)
nClustGap = maxSE(gap_statKmeans$Tab[,3], 
                  gap_statKmeans$Tab[,4], method = "firstSEmax")


#Elbow method:
clusterVec = 1:20 # must begin with 1
logLikeVec = rep(NA, length(clusterVec))
wcsseVec = rep(NA, length(clusterVec))
for (k in clusterVec) { 
  kfit = kmeans(projectionNorm, centers = k)
  logLikeVec[k] = getLogLike(assay(rse), kfit$cluster)
  wcsseVec[k] = kfit$tot.withinss
}
logLikeElbow = getElbow(logLikeVec, clusterVec)
wcsseElbow = getElbow(wcsseVec, clusterVec)

#repeat for WCCSE:


# elbowMethod:

getElbow = function(measureVec, clusterVec){
  elbowMetric  = rep(NA, length(clusterVec))
  for (k in clusterVec) {
    fit1 = lm(measureVec[1:k] ~ clusterVec[1:k])
    fit2 = lm(measureVec[k:length(clusterVec)] ~ clusterVec[k:length(clusterVec)])
    elbowMetric[k] = sum(fit1$residuals^2) + sum(fit2$residuals^2)
  }
  return(elbowMetric)
}

out = data.frame(nPCs = PCrange, 
  type = 
  c("Distortion", "Silhouette", "GapStat", "logLikeElbow", "wcsseElbow"),
           nClusters = c(
              which.min(nClustDistortion), 
              which.max(nClustSilhouette), 
              nClustGap,
              which.min(logLikeElbow),
              which.min(wcsseElbow))
)
write.csv(out, file.path(peaksDir, paste0(sampleName, PCrange, "NClusters.csv")),
          row.names = F)


