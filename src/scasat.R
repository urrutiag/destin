#!/usr/bin/env Rscript

#SBATCH --mem=32g
#SBATCH --time=24:00:00

args = commandArgs(trailingOnly=TRUE)

sampleName = args[1]
outDir = args[2]
cellDataFile = args[3]

peaksDir = file.path(outDir, "peaks")

library(SummarizedExperiment)
library(proxy)
library(cluster)

#load data:
fpath <- file.path(peaksDir, paste0(sampleName, "500KPeaksAnnotated.Rdata"))
load(fpath)

nClusters = length(unique(rse$cell_type))

### keep top nPeaks regions

qValCutoff = -log10(.2)
peakBelowCutoff = rowRanges(rse)$qValue >= qValCutoff # qvalue on -log scale
if (!is.null( qValCutoff))  nPeaks = sum(peakBelowCutoff)
rse = rse[peakBelowCutoff,]

### QC ------------------------------------------------------------
  
rse = rse[Matrix::rowSums(assay(rse)) >= 10, ]
cellSum = Matrix::colSums(assay(rse))
cutoffs = 2 ^ ( median(log2(cellSum)) + c(-3,3)*mad(log2(cellSum)) )
rse = rse[ , cellSum > cutoffs[1] & cellSum < cutoffs[2] ]
 
jacDist = proxy::dist(as(t(assay(rse)), "matrix"), method = "Jaccard")
kfitScasat = cluster::pam(jacDist, k = nClusters, diss=TRUE)

purity = sum( apply(table(kfitScasat$cluster, rse$cell_type), 1, max) ) / ncol(rse)

summary = data.frame(method = "scasat",
                     sampleName = sampleName,
                      nPeaksActual = nrow(rse),
                      purity = purity
          )

write.csv(summary,
 file.path(peaksDir, paste0(sampleName, "ScasatResults.csv")),
 row.names = F)

