\name{plotNClusters}
\alias{plotNClusters}
\title{
Create a plot of estimated number of clusters.
}
\description{
Plot the results of estimateNClusters().  Estimate is shown in blue by method.
}
\usage{
plotNClusters(clusterEst)
}

\arguments{
  \item{clusterEst}{ results of estimateNClusters()}
}

\value{
plot is printed 
}

\author{
 Gene Urrutia \email{urrutia@email.unc.edu} 

}

\examples{

nClustersMax = 20
clusterEst = estimateNClusters(rse, nClustersMax = nClustersMax)

plotNClusters(clusterEst)
                      

\keyword{ package }