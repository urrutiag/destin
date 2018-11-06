\name{plotClusterTsne}
\alias{plotClusterTsne}

\title{
Create t-SNE plot of cluster results
}
\description{
Performs t-SNE clustering and plots. For visualization only.  T-SNE is not used for Destin clustering. 

Input are weighted principle components after tuning parameters to optimized to maximize mulitnomial likelihood.  Cells are colored according to Destin clustering.  Can optionally label the clusters. 
}
\usage{
plotClusterTsne(clusterResults, clusterLabels = NULL)
}

\arguments{
  \item{clusterResults}{ results of destinGrid()}
  \item{clusterLabels}{ optional vector of cluster labels }
}

\value{
plot is printed 
}

\author{
 Gene Urrutia \email{urrutia@email.unc.edu} 

}

\examples{

results = destinGrid (rse, sampleName,
                      PCrange = PCrange,
                      TSSWeightsList = TSSWeightsList,
                      DHSWeightsList = DHSWeightsList,
                      nClusters = nClusters,
                      nCores = nCores)
                      
plotClusterTsne(results)

\keyword{ package }
