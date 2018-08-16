\name{estimateNClusters}
\alias{estimateNClusters}
\title{
Estimate the number of clusters
}
\description{
Models the cluster fit using multinomial likelihood and selects optimal number of cluster based on systematic elbow method
}
\usage{
estimateNClusters(rse, clusterVec = 1:20)
}
\arguments{
    \item{ rse }{
    ranged summarized experiment containing region by cell accessiblity matrix, region data, and cell data}
    \item{ clusterVec }{
    range of number of clusters to test }
}

\value{
list:

logLikes:  data.frame number of clusters and model based log likelihood 

nClustersEstimate: number of clusters pertaining to the best fit
}

\author{
 Gene Urrutia \email{urrutia@email.unc.edu} 
}
\examples{
library(destin)

set.seed(10)
clusterEst = estimateNClusters(rse,   clusterVec = 1:5)
qplot(data = clusterEst$logLikes, x = nClusters, y = logLike) 
clusterEst$nClustersEstimate
}
\keyword{ package }
