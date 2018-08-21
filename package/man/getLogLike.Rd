\name{getLogLike}
\alias{getLogLike}
\title{
  Calculate the log likelihood of a cluster result.
}
\description{
  Calculates the log likelihood of a cluster result based on the multinomial model.  
}
\usage{
getLogLike(countMat, cluster)
}
\arguments{
  \item{ countMat }{
    accessibility matrix  (region by cell)
  }
  \item{ cluster }{
    vector of cluster membership matching the cells of the accessibility matrix 
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
  
  logLike =  getLogLike(countMat, kfit$cluster)
  
}
\keyword{ package }
