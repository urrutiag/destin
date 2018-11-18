\name{getGeneAccessibility}
\alias{getGeneAccessibility}
\title{
calculates the proportion of total accessibility by cell type
}
\description{
This step calculates the proportion of total accessibility by cell type.  This is required input for magma and EWCE.

The result is a stacked data frame containing the cell type specificity for each gene - cell type combination.  cell type specificity is calculated as the ratio between the cell type accessibility and the total accessibility across all cell types.
}
\usage{
getGeneAccessibility(se)
}
\arguments{
  \item{se}{
 summarized experiment of chromatin accessiblility aggregated by gene
}
}
\value{
data frame containing the cell type specificity for each gene - cell type combination.  cell type specificity is calculated as the ratio between the cell type accessibility and the total accessibility across all cell types.
}
\author{
 Gene Urrutia \email{gene.urrutia@gmail.com} 
}
\examples{
geneAccessibility = getGeneAccessibility(se)
}
\keyword{ package }