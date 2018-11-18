\name{getQuantsForMagma}
\alias{getQuantsForMagma}
\title{
convert the cell type specific accessibility from the geneAccessibility data frame to 40 binned quantiles
}
\description{
converts the cell type specific accessibility from the geneAccessibility data frame to 40 binned quantiles

As done by Skeene et al., we binned the cell type specificity into 40 quantiles, with an additional bin for regions with no accessibility. The cell type specific bin was used as the feature vector for association testing. Finally test for trend between cell type specificity was performed. 

The specificity quantile matrix must be printed to file, as it is required by magma using the terminal (shell) command line.
}

\usage{
getQuantsForMagma( geneAccessibility)
}
\arguments{
  \item{geneAccessibility}{
the cell type specific accessibility data frame created by getGeneAccessibility
}
}
\value{
data frame one row per gene, and one column for each cluster, with a value from 0 to 40.  Genes in which a cluster has high cellTypeSpecificity receives a high quantile, and vice versa.  The highest quantile is 40 and 0 indicates that there was no accessibility for any cell in the cluster.
}
\author{
 Gene Urrutia \email{gene.urrutia@gmail.com} 
}
\examples{
quantsForMagma = getQuantsForMagma( geneAccessibility)
}
\keyword{ package }
