\name{getDestin}
\alias{getDestin}
\title{
Cluster individual cells by chromatin accessibility using region weighted PCA and k-means clustering.
}
\description{
Cluster individual cells by chromatin accessibility using region weighted PCA and k-means clustering..  The function takes a single choice of DHS weights and THS weights,  as well as a range of number of principle components.  
}
\usage{
getDestin( rse, nClusters, PCrange=10, TSSWeights=c(1,1), DHSWeights=c(1,1)) }
\arguments{
    \item{ rse }{
    ranged summarized experiment containing region by cell accessiblity matrix, region data, and cell data}
    \item{ PCrange }{
    range of principle components to test }
    \item{ nClusters }{
    nClusters as calculated by \code{\link{estimateNClusters} } }
    \item{ TSSWeights }{
    Single choice of transcription start site weights. First weight for promoter region, second weight for distal element}
    \item{ DHSWeights }{
    Single choice of DNase hypersensitivity site weights. Weights correspond to weights 1 and 2 of beta distribution, mapping database accessibility to a region weight.}
    \item{ ... }{
      parameters for development use
    }
}

\value{
Results of destin are a list for each nPC in selected range.  Each list contains:

Summary: data.frame containing one row, the result of the optimal choice of hyperparameters based on the multinomial likelihood.  The row contains DHS weight choice, THS weight choice, number of principle components, number of regions, and model based loglikelihood. 

Cluster: the clustering result.

PCs: Principle component projections (cell x nPCs) where the number of PCs is determined by maximum likelihood
}

\author{
 Gene Urrutia \email{gene.urrutia@gmail.com} 
}
\examples{

results = getDestin( rse, nClusters = nClusters)

}
\keyword{ package }
