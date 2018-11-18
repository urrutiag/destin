\name{getEWCE}
\alias{getEWCE}
\title{
determine whether GWAS results are associated with increased accessibility in a particular cell type cluster
}
\description{
Expression Weighted Celltype Enrichment (ECWC) (Skene and Grant, 2016).  ECWC tests whether the set of significant SNPs are more specific to cell type than are a random set of SNPs. 

The empirical score = sum of cellTypeSpecificity from the geneAccessibility data frame only within the snpGenes, specific by cell type.

The p-value is calculated via permutation across all genes of the above score.

getEWCE is a function which takes as input: 
- geneAccessibility: the cell type specific accessibility data frame created is the previous step
- geneList:  list of significant snps for a single GWAS
}
\usage{
getEWCE( geneAccessibility, geneList, nBoot = 1000 )
}
\arguments{
  \item{geneAccessibility}{
the cell type specific accessibility data frame created by getGeneAccessibility
}
  \item{geneList}{
list of significant snps from a GWAS
}
}
\value{
data frame with p-value, foldChange, sd from mean for each cluster (cell type)
}
\author{
 Gene Urrutia \email{gene.urrutia@gmail.com} 
}
\examples{
resultsEWCE = getEWCE(geneAccessibility, geneList)
}
\keyword{ package }
