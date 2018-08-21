\name{getDiffAccess}
\alias{getDiffAccess}
\title{
Calculate cluster-specific differentially accessible regions 
}
\description{
Calculates cluster-specific differentially accessible regions by Fisher's exact test. Parallel cpu computing is supported.  Each cluster is tested against all other clusters combined.
}
\usage{
getDiffAccess(rse, clusterMap, nCores = NULL)
}
\arguments{
    \item{ rse }{
    ranged summarized experiment containing region by cell accessiblity matrix, region data, and cell data}
    \item{ clusterMap }{
    data.frame of cluster identity (cellID and cluster) as output by 
    \code{\link{destinGrid} }  }
   \item{ nCores }{
    number of CPU cores to use for parallel computing } 
}

\value{
Results are appended to the row ranges of ranged summarized experiment and is meant to replace input. Specifically, results include p-value, Bonferroni corrected p-value, and log2 fold change (+1 pseudo-count).  
}

\author{
 Gene Urrutia \email{urrutia@email.unc.edu} 
}
\examples{
library(destin)

rse = getDiffAccess(rse, clusterMap = results$cluster)

}
\keyword{ package }
