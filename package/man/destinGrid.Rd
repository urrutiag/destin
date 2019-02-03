\name{destinGrid}
\alias{destinGrid}
\title{
  Cluster individual cells by chromatin accessibility using region weighted PCA and k-means clustering.
}
\description{
  Cluster individual cells by chromatin accessibility using region weighted PCA and k-means clustering..  The function takes multiple choices of DHS weights and THS weights, as well as a range of number of principle components.  These hyperparameters are optimized by multinomial log likelihood and the optimal clustering result is returned.  Parallel cpu computing is supported.
}
\usage{
  getDestin(rse, nClusters, 
            PCrange = 3:25,
            TSSWeightsList = list(c(1,2), c(1,1.5), c(1,1.25), c(1,1)),
            DHSWeightsList = list(c(1,1), c(1,2), c(1,3), c(1,5)),
            nCores=NULL, ...) }
\arguments{
  \item{ rse }{
    ranged summarized experiment containing region by cell accessiblity matrix, region data, and cell data}
  \item{ nClusters }{
    nClusters as calculated by \code{\link{estimateNClusters} } }
  \item{ PCrange }{
    range of principle components to test }
  \item{ TSSWeights }{
    List of mulitple transcription start site weights. First weight for promoter region, second weight for distal element}
  \item{ DHSWeights }{
    List of mulitple DNase hypersensitivity site weights. Weights correspond to weights 1 and 2 of beta distribution, mapping database accessibility to a region weight.}
  \item{ nCores }{
    number of CPU cores to use for parallel computing } 
  \item{ ... }{
    parameters for development use
  }
}

\value{
Results of destinGrid are a list:

Summary: data.frame containing one row, the result of the optimal choice of hyperparameters based on the multinomial likelihood.  The row contains DHS weight choice, THS weight choice, number of principle components, number of regions, and model based loglikelihood. 

Cluster: the clustering result.

PCs: Principle component projections (cell x nPCs) where the number of PCs is determined by maximum likelihood

SummaryAll: As summary above, except containing results for every weight / nPC combination.  

}

\author{
  Gene Urrutia \email{gene.urrutia@gmail.com} 
}
\examples{
  
results = destinGrid (rse, nClusters = nClusters)

}
\keyword{ package }
