\name{plotClusterTsne}
\alias{plotClusterTsne}

\title{
Create PCA or t-SNE plot of cluster results
}
\description{
Either plots PCA or performs t-SNE clustering and plots.  T-SNE is highly variable, not used for Destin clustering, and is for visualization only. 

Input are weighted principle components after tuning parameters to optimized to maximize mulitnomial likelihood.  Cells are colored according to Destin clustering.  Can optionally label the clusters. 
}
\usage{
plotCluster(clusterResults, type = "PCA", clusterLabels = NULL)
}

\arguments{
  \item{clusterResults}{ results of destinGrid()}
  \item{clusterLabels}{ optional vector of cluster labels }
}

\value{
plot is printed 
}

\author{
 Gene Urrutia \email{gene.urrutia@gmail.com} 

}

\examples{

results = destinGrid (rse, nClusters = nClusters,
                      nCores = nCores)
                      
plotCluster(results)

}
\keyword{ package }
