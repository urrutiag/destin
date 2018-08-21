\name{destinGrid}
\alias{destinGrid}
\title{
  Cluster individual cells by chromatin accessibility using region weighted PCA and k-means clustering.
}
\description{
  Cluster individual cells by chromatin accessibility using region weighted PCA and k-means clustering..  The function takes multiple choices of DHS weights and THS weights, as well as a range of number of principle components.  These hyperparameters are optimized by multinomial log likelihood and the optimal clustering result is returned.  Parallel cpu computing is supported.
}
\usage{
  getDestin(rse, sampleName, 
            PCrange = 3:25,
            TSSWeightsList = list(c(1,2), c(1,1.5), c(1,1.25), c(1,1)),
            DHSWeightsList = list(c(1,1), c(1,2), c(1,3), c(1,5)),
            nClusters, nCores=NULL, writeOut=F, outDir=NULL) }
\arguments{
  \item{ rse }{
    ranged summarized experiment containing region by cell accessiblity matrix, region data, and cell data}
  \item{ PCrange }{
    range of principle components to test }
  \item{ TSSWeights }{
    List of mulitple transcription start site weights. First weight for promoter region, second weight for distal element}
  \item{ DHSWeights }{
    List of mulitple DNase hypersensitivity site weights. Weights correspond to weights 1 and 2 of beta distribution, mapping database accessibility to a region weight.}
  \item{ nClusters }{
    nClusters as calculated by \code{\link{estimateNClusters} } }
  \item{ nCores }{
    number of CPU cores to use for parallel computing } 
  \item{ writeOut }{
    option to write results to disk } 
  \item{ outDir }{
    if writeOut is set to TRUE, this is the destination directory, current working directory by default } 
}

\value{
Results of destin are a list:

The first item is a data.frame containing one row, the result of the optimal choice of hyperparameters based on the multinomial likelihood.  The row contains DHS weight choice, THS weight choice, number of principle components, number of regions, purity if cell type is known, and model based loglikelihood. 

The second item is the clustering result. 
}

\author{
  Gene Urrutia \email{urrutia@email.unc.edu} 
}
\examples{
library(destin)
  
clusterEst = estimateNClusters(rse, clusterVec = 1:5)
nClusters = clusterEst$nClustersEstimate
PCrange = 3:10
TSSWeightsList = list(c(1,2), c(1,1.5), c(1,1.25), c(1,1))
DHSWeightsList = list(c(1,1), c(1,2), c(1,3), c(1,5))
results = destinGrid (rse, sampleName, 
            PCrange = PCrange,
            TSSWeightsList = TSSWeightsList,
            DHSWeightsList = DHSWeightsList,
            nClusters = nClusters)
results
}
\keyword{ package }
