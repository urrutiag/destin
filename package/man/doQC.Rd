\name{doQC}
\alias{doQC}
\title{
Remove cells and regions that do not pass quality control
}
\description{
Removes cells and regions that do not pass qualtiy control.  Specifically, removes regions not detected by at least 5 cells.  Removed cells that do not fall within 3 standard deviations of median (calculated robustly using median absolute deviations).
}
\usage{
doQC(rse, regionSumCutoff = 5, cellSumCutoffSDs = 3)
}
\arguments{
    \item{ rse }{
    ranged summarized experiment containing region by cell accessiblity matrix, region data, and cell data }
    \item{ regionSumCutoff }{
    exclude regions whose total accessibility count is less than regionSumCutoff
}
    \item{ cellSumCutoffSDs }{
    exclude cells whose total accessibility count is outside of cellSumCutoffSDs number of sds away from median 
}

}

\value{
Result is a QC'd ranged summarized experiment and is meant to replace input.   

}

\author{
 Gene Urrutia \email{gene.urrutia@gmail.com} 
}
\examples{
library(destin)

rse = doQC(rse, regionSumCutoff = 5, cellSumCutoffSDs = 3)


}
\keyword{ package }
