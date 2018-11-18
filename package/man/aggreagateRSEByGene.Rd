\name{aggregateRSEByGene}
\alias{aggregateRSEByGene}
\title{
Aggregate chromatin accessible peaks by gene.
}
\description{
The ranged summarized experiment contains multiple accessible peaks per entrez ID.  Thus we aggregate across peaks within gene to create a chromatin accessiblity matrix that is no longer binary.

This step is high processing, thus it is recommended to run in parallel, but that is optional.  It is also high memory, thus we use sparseMatrix format.

The output is an aggregated summarized experiment where rows are genes by entrez ID
}
\usage{
aggregateRSEByGene(rse, nCores = NULL)
}
\arguments{
  \item{rse}{
ranged summarized experiment
}
  \item{nCores}{
cpu cores to use in parallel
}
}
\value{
summarized experiment where rows are genes by entrez ID
}
\author{
 Gene Urrutia \email{gene.urrutia@gmail.com} 
}
\examples{
se = aggregateRSEByGene(rse)
}
\keyword{ package }