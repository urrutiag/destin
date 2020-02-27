\name{annotateHuman}
\alias{annotateHuman}
\title{
Annotate human Ensemble features from a ranged summarized experiment with human Entrez ID and HGNC symbol.
}
\description{
Utilizes \href{https://www.rdocumentation.org/packages/ChIPpeakAnno/versions/3.6.5/topics/addGeneIDs}{addGeneIDs} function from \href{https://www.rdocumentation.org/packages/ChIPpeakAnno/versions/3.6.5}{ChIPpeakAnno} package to annotate human Entrez ID and HGNC symbol from Ensemble feature. \href{https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html}{org.Hs.eg.db data} set is used for annotation. 

annotated Esememble peak feature (human) -> HGNC symbol (human) + Entrez ID (human) 
}
\usage{
annotateHuman(rse)
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
library(destin)

# build toy rse to test functionality
rse = genTestRSE()

# adds Ensembl feature using ChIPpeakAnno's annotatePeakInBatch().  This is part of the standard Destin workflow.
rse = annotateRSE(rse, "hg19")
  
# annotate entrez id using using ChIPpeakAnno's addGeneIDs().  
rse = annotateHuman(rse)
  
# aggregate (modified to accomodate null mgi_symbol)

# Test downstream workflow
# Aggregate chromatin accessible peaks by gene.
se = aggregateRSEByGene(rse)
# calculates the proportion of total accessibility by cell type
getGeneAccessibility(se)

}
\keyword{ package }
