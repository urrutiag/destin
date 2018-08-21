\name{doQC}
\alias{doQC}
\title{
Remove cells and regions that do not pass quality control
}
\description{
Removes cells and regions that do not pass qualtiy control.  Specifically, removes regions not detected by at least 5 cells.  Removed cells that do not fall within 3 standard deviations of median (calculated robustly using median absolute deviations).
}
\usage{
doQC(rse)
}
\arguments{
    \item{ rse }{
    ranged summarized experiment containing region by cell accessiblity matrix, region data, and cell data }
}

\value{
Result is a QC's ranged summarized experiment and is meant to replace input.   

}

\author{
 Gene Urrutia \email{urrutia@email.unc.edu} 
}
\examples{
library(destin)

rse = doQC(rse)


}
\keyword{ package }
