\name{getLogLike}
\alias{getLogLike}
\title{
  Calculate the log likelihood of a cluster result.
}
\description{
  Calculates the log likelihood of a cluster result based on the multinomial model.  
}
\usage{
getLogLike(countMat, cluster, sum = T)
}
\arguments{
  \item{ countMat }{
    accessibility matrix  (region by cell)
  }
  \item{ cluster }{
    vector of cluster membership matching the cells of the accessibility matrix 
  }
  \item{ sum }{
    if TRUE, output is the sum of log likelihoods across all cells according to assigned cluster membership. 
    if FALSE, output is a cell by cluster matrix, with each entry the log likelohood of the cell given cluster.
  }  
}
\value{
  Log likelihood of a clustering result based on multinomial model
}

\author{
  Gene Urrutia \email{urrutia@email.unc.edu} 
}
\examples{
  library(destin)
  
  # returns the loglikelihood 
  logLike =  getLogLike(countMat, kfit$cluster)
  
}
\keyword{ package }
