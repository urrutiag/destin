\name{createRSEfrom10xMatrix}
\alias{createRSEfrom10xMatrix}
\title{
Create ranged summarized experiment from 10x peak by cell matrix
}
\description{
Create ranged summarized experiment from 10x peak by cell matrix
}
\usage{
rse = createRSEfrom10xMatrix(data10xDir)
}
\arguments{
    \item{ data10xDir }{
    path to directory which contains 10x peak by cell matrix files: 
    1) barcodes.tsv
    2) matrix.mtx
    3) peaks.bed
    }
}

\value{
ranged summarized experiment containing region by cell accessiblity matrix, region data, and cell data.
}

\author{
 Gene Urrutia \email{gene.urrutia@gmail.com} 
}
\examples{
library(destin)

rse = createRSEfrom10xMatrix(data10xDir)

}
\keyword{ package }

