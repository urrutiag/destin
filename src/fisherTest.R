#!/usr/bin/env Rscript

#SBATCH mem=32g
#SBATCH time=24:00:00
#SBATCH cpus-per-task=24

args = commandArgs(trailingOnly=TRUE)

sampleName = args[1]
outDir = args[2]

# sampleName = "BuenrostroMouse"
# outDir = "/proj/yuchaojlab/gene/scATACseq/data/BuenrostroMouse"

library(SummarizedExperiment)
library(data.table)
library(parallel)

peaksDir = file.path(outDir, "peaks")

clusterFile = dir(peaksDir, pattern="FinalCluster", full.names = T)
dataFile = dir(peaksDir, pattern="500KPeaks", full.names = T)

load(dataFile)

### QC ------------------------------------------------------------
rse = rse[Matrix::rowSums(assay(rse)) >= 5, ]
cellSum = Matrix::colSums(assay(rse))
cutoffs = 2 ^ ( median(log2(cellSum)) + c(-3,3)*mad(log2(cellSum)) )
rse = rse[ , cellSum > cutoffs[1] & cellSum < cutoffs[2] ]

clusterData = fread(clusterFile)

# all.equal(clusterData$cellID, rse$cellID)

#INPUT:
C = clusterData$cluster
X = assay(rse)

#table(C, rse$cell_type)
#C[C==5] = 1

# 4 sec / 1000 loci 

nCores = 24
cl = makeCluster(nCores)
clusterEvalQ(cl, library(Matrix))
clusterExport(cl, list("C", "X")) 

pfcList = list()
for (cellIndex in seq_along(unique(C))) {

Cref = C == unique(C)[cellIndex]
freqRef = ( Matrix::rowSums(X[,Cref]) + 1) / sum(Cref)
freqRest = ( Matrix::rowSums(X[,!Cref]) + 1) / sum(!Cref)
log2FC = log2(freqRef / freqRest)

clusterExport(cl, list("Cref")) 
pVals = parSapply(cl, 1:nrow(X), function(regionIndex) {
  mat = table(Cref, X[regionIndex,])
  fish = fisher.test(mat)
  return( fish$p.value )
} )

pValsCor = p.adjust(pVals, method = "bonferroni")

out = data.frame(pVals,
          pValsCor,
          log2FC)

names(out) = paste( c("pVals", "pValsCor", "log2FC"), 
   unique(C)[cellIndex], sep = "_")

pfcList[[cellIndex]] = out

}

stopCluster(cl)

pfc = do.call("cbind", pfcList)
pfc = cbind(pfc)
colData(rse)$cluster = C
elementMetadata(rse) = cbind(elementMetadata(rse), pfc)

write.csv(rowRanges(rse), 
  file = file.path(peaksDir, paste0(sampleName, "Fisher.csv")),
  row.names = F)

