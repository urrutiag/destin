\name{annotateMouseToHuman}
\alias{annotateMouseToHuman}
\title{
Annotate genes from a ranged summarized experiment with human Entrez ID and HGNC symbol via MGI symbol starting with mouse annotated peak feature .
}
\description{
GWAS predominantly use human subjects, while our chromatin acessibility data set is mouse forebrain.  Thus we need to utilize a series of 3 conversions to map the accessible regions to HGNC and Entrez ID.

annotated peak feature (mouse) -> MGI symbol  (mouse) -> HGNC symbol (human) -> Entrez ID (human) 
}
\usage{
annotateMouseToHuman(rse)
}
\arguments{
  \item{rse}{
ranged summarized experiment
}
}
\value{
ranged summarized experiment with human gene annotations
}
\author{
 Gene Urrutia \email{gene.urrutia@gmail.com} 
}
\examples{
rse = annotateMouseToHuman(rse)
}
\keyword{ package }
