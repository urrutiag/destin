\name{getElbow}
\alias{getElbow}
\title{
Systematically determine elbow of plot
}
\description{
Fits a linear spline with a single knot to each point on a plot and return the choice with loest sum of squared errors. This is used to determine optimal number of clusters.
}
\usage{
getElbow(measureVec, clusterVec)
}
\arguments{
    \item{ measureVec }{
    measure (log likelihood) to model as a function of number of clusters }
    \item{ clusterVec }{
    range of number of clusters to test }
}

\value{
Vector of sum squared errors for each value of knot.  The lowest value is selected as optimal and estimates the optimal number of clusters.
}

\author{
 Gene Urrutia \email{urrutia@email.unc.edu} 
}
\examples{
library(destin)

  logLikeVec = rep(NA, length(clusterVec))
  for (k in clusterVec) { 
    kfit = kmeans(projectionNorm, centers = k)
    logLikeVec[k] = getLogLike(assay(rse), kfit$cluster)
  }
  
  logLikeElbow = getElbow(logLikeVec, clusterVec)
  
  nClustersEstimate =  which.min(logLikeElbow)


}
\keyword{ package }
