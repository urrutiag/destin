\name{getLogLike}
\alias{getLogLike}
\title{
  Calculate the log likelihood of a cluster result.
}
\description{
  Calculates the log likelihood of a cluster result based on the multinomial model.  
}
\usage{
getLogLike(countMat, cluster, reassign = F)
}
\arguments{
  \item{ countMat }{
    accessibility matrix  (region by cell)
  }
  \item{ cluster }{
    vector of cluster membership matching the cells of the accessibility matrix 
  }
  \item{ reassign }{
    if FALSE, output is the sum of log likelihoods across all cells according to assigned cluster membership. 
    if TRUE, output is a cell by cluster matrix, with each entry the log likelohood of the cell belonging to that cluster.  
  }  
}
\value{
  Log likelihood of a clustering result based on multinomial model
}

\author{
  Gene Urrutia \email{gene.urrutia@gmail.com} 
}
\examples{
  library(destin)
  
  # returns the loglikelihood 
  logLike =  getLogLike(countMat, kfit$cluster)
  
}
\keyword{ package }
