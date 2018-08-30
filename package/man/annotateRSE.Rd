\name{annotateRSE}
\alias{annotateRSE}
\title{
Add region annotation
}
\description{
Adds region annotation to a ranged summarized experiment. Annotation includes distance to transcription start site, annotation based on the distance, and mean accessibility accross a refererence set of DNase hypersensitivity data. 
}
\usage{
annotateRSE(rse, model)
}
\arguments{
    \item{ rse }{
    ranged summarized experiment containing region by cell accessiblity matrix, region data, and cell data}
    \item{ model }{
    model used for read alignment. hg19 or mm10 supported }
}

\value{
Ranged summarized experiment containing row ranges with annotated regions.  Annotations: DHSsum (sum accessibilty of reference set), distanceToTSS, feature, region type (promoter or distal element)
}

\author{
 Gene Urrutia \email{urrutia@email.unc.edu} 
}
\examples{
library(destin)

rse = annotateRSE(rse, model)

}
\keyword{ package }
