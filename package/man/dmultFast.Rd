\name{dmultFast}
\alias{dmultFast}
\title{
Calculates the log likelihood of a cluster result for an individual cell.
}
\description{
Calculates the log likelihood of a cluster result for an individual cell.  The function is optimized for binary accessibility.  
}
\usage{
dmultFast(x, prob)
}
\arguments{
    \item{ x }{
    accessibility vector for individual cell
}
    \item{ prob }{
    mean accessibility of all cells belonging to same cluster as cell}
}

\value{
Log likelihood of individual cell based on multinomial model
}

\author{
 Gene Urrutia \email{gene.urrutia@gmail.com} 
}
\examples{
library(destin)

logLike = dmultFast(x = countMat[,myCellIndex],
              prob = empiricalProbList[[paste(cluster[myCellIndex])]])

}
\keyword{ package }
