installed <- rownames(installed.packages())
pkgs = c("cluster", "data.table", "ggplot2",
         "gridExtra", "irlba",  "Matrix",
         "parallel", "Rtsne")
pkgs <- setdiff(pkgs, installed)
if (length(pkgs))
  install.packages(pkgs, dep=c("Depends", "Imports"))

biocPkgs = c("ChIPpeakAnno", "GenomicAlignments", "rtracklayer")
biocPkgs <- setdiff(biocPkgs, installed)
if (length(biocPkgs)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(biocPkgs)
}

#ClusterR is an optional package:
if ( ! "ClusterR" %in% rownames(installed.packages() ) )
  install.packages("ClusterR", dep=c("Depends", "Imports"))


yourPathToDestinRepo = "~/Documents/gitRepos/destin"
install.packages(file.path(yourPathToDestinRepo,"package"), repos = NULL, type = "source")

install.packages("proxy")


install.packages("devtools")
library(devtools)
devtools::install_github("SUwonglab/scABC")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("BiocManager")
BiocManager::install("chromVAR")
biocPkgs = c("motifmatchr", "BSgenome.Hsapiens.UCSC.hg19",
             "BSgenome.Mmusculus.UCSC.mm10", "JASPAR2016")
BiocManager::install(biocPkgs)

