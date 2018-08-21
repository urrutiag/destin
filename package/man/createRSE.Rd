\name{createRSE}
\alias{createRSE}
\title{
Create ranged summarized experiment from peak file and bam files
}
\description{
Create ranged summarized experiment from peak file and bam files.
}
\usage{
rse = createRSE(bamDir, bamFiles, bedData)
}
\arguments{
    \item{ bamDir }{
    path to directory which contains sorted and indexed bam files
    }
    \item{ bamFiles }{
    vector of file names of sorted and indexed bam files
    }
    \item{ bedData }{
    peaks data in GRanges format, see example for reading MACS2 narrow peaks file
    }
}

\value{
ranged summarized experiment containing region by cell accessiblity matrix, region data, and cell data.
}

\author{
 Gene Urrutia \email{urrutia@email.unc.edu} 
}
\examples{
library(destin)

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
bedData = rtracklayer::import(bedFile, format = "BED",
               extraCols = extraCols_narrowPeak)

bamFiles = dir( bamDir, pattern = "final.bam" )
bamFiles = bamFiles[!grepl("bai", bamFiles)]

rse = createRSE(bamDir, bamFiles, bedData)

}
\keyword{ package }
