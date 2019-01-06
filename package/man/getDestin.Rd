\name{getDestin}
\alias{getDestin}
\title{
Cluster individual cells by chromatin accessibility using region weighted PCA and k-means clustering.
}
\description{
Cluster individual cells by chromatin accessibility using region weighted PCA and k-means clustering..  The function takes a single choice of DHS weights and THS weights,  as well as a range of number of principle components.  
}
\usage{
getDestin( rse, PCrange=10, TSSWeights=c(1,1), DHSWeights=c(1,1), 
                     nClusters, outCluster = F ) }
\arguments{
    \item{ rse }{
    ranged summarized experiment containing region by cell accessiblity matrix, region data, and cell data}
    \item{ PCrange }{
    range of principle components to test }
    \item{ TSSWeights }{
    Single choice of transcription start site weights. First weight for promoter region, second weight for distal element}
    \item{ DHSWeights }{
    Single choice of DNase hypersensitivity site weights. Weights correspond to weights 1 and 2 of beta distribution, mapping database accessibility to a region weight.}
    \item{ nClusters }{
    nClusters as calculated by \code{\link{estimateNClusters} } }
}

\value{
Results of destin are a data.frame containing one row for each choice of number of principle components.  Each row contains DHS weight choice, THS weight choice, number of principle components, number of regions, purity if cell type is known, and model based loglikelihood.  If outCluster is set to TRUE, then the output is a list where the first item is the above, and the second item is the clustering result. 
}

\author{
 Gene Urrutia \email{gene.urrutia@gmail.com} 
}
\examples{
library(destin)

clusterEst = estimateNClusters(rse, clusterVec = 1:5)
nClusters = clusterEst$nClustersEstimate
results = getDestin( rse, nClusters = nClusters, outCluster = T)
results

}
\keyword{ package }
