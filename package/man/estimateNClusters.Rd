\name{estimateNClusters}
\alias{estimateNClusters}
\title{
Estimate the number of clusters
}
\description{
Destin uses an unsupervised model-based likelihood to estimate number of clusters. Like WCSSE, the likelihood is non-decreasing with increasing number of clusters.  We fit the model log-likelihood as a linear spline function of number of clusters with a single knot. The knot placement which results in the best fit (lowest sum of squared errors) is selected as number of clusters (n_k).

Also calculated using 4 other methods to estimate number of clusters:

wcsseElbow: They fit WCSSE as a linear spline function of k with a single knot. The knot placement which results in the best fit (lowest sum of squared errors) is selected as n_k

silhouette: n_k is selected as the number of clusters that maximizes average silhouette width across all cells

distortion: accounts for the decrease in WCSSE by comparing to a baseline model incorporating a weight factor based on number of dimensions

GapStat: The gap (Tibshirani et al., 2001) statistic accounts for the tendency to overestimate n_k by using a permutation approach and fitting a null distribution to the WCSSE which naturally decreases as n_k increases. Then the empirical WCSSE is compared to the null distribution with the gap statistic as output. To select the number of clusters the “first SE” method is used which determines the first nk where the gap statistic is within 1 standard deviation of the following gap statistic.
}
\usage{
clusterEst = estimateNClusters(rse, nClustersMax)
}
\arguments{
    \item{ rse }{
    ranged summarized experiment containing region by cell accessiblity matrix, region data, and cell data}
    \item{ nClustersMax }{
    maximum number of clusters to test, search space begins at 1 }
}

\value{
list:

logLikes:  data.frame number of clusters and model based log likelihood 

nClustersEstimate: number of clusters pertaining to the best fit
}

\author{
 Gene Urrutia \email{urrutia@email.unc.edu} 
}
\examples{
library(destin)

nClustersMax = 20
clusterEst = estimateNClusters(rse, nClustersMax = nClustersMax)

plotNClusters(clusterEst)
}
\keyword{ package }
