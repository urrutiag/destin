
yourPathToDestinRepo = "~/Documents/gitRepos/destin"
install.packages(file.path(yourPathToDestinRepo,"package"), repos = NULL, type = "source")
library(destin)
library(gridExtra)

dataDir = '~/Dropbox/Gene_shared/scATACseq/destin/data'
sampleNames = c('BuenrostroHuman', 'BuenrostroMouse', 'Corces500',
                'CorcesSub', 'CusanovichGM12878vsHEK', 'CusanovichGM12878vsHL',
                'CusanovichGM12878vsPatski', 'PreisslP56Annotated')
models = c('hg19', 'mm10', 'hg19', 'hg19', 'hg19', 'hg19', 'hg19', 'mm10')
nClustersList = c(6,2,4,13,2,2,2,8)
cellTypeNames = c('cell_type', 'cell_type', 'cell_type', 'cell_type', 'cellType', 'cellType', 'cellType', 'cell_type')

# cluster purity

for (dataIndex in 1:8) {
  getPurityByDepthAdjustment(dataIndex = dataIndex)
}


getPurityByDepthAdjustment = function(dataIndex){
  
  nClusters = nClustersList[dataIndex]
  sampleName = sampleNames[dataIndex]
  model = models[dataIndex]
  cellTypeName = cellTypeNames[dataIndex]
  
  filePath = dir(dataDir, full.names = TRUE, pattern = sampleName)
  load(filePath) 
  assay(rse) = as(assay(rse), "dgCMatrix")
  rse$cellID = rse[[cellTypeName]]
  rse = annotateRSE(rse, model)
  rse = doQC(rse, regionSumCutoff = 5, cellSumCutoffSDs = 3)
  
  depthAdjustments = c("postPCA", "prePCA", "none")
  pcaTableList = lapply( depthAdjustments, function(depthAdjustment) {
    results = destinGrid (rse, nClusters, nCores = 7, depthAdjustment = depthAdjustment)
    clusterTable = table(results$cluster$cluster, results$cluster$cellID)
    clusterPurity = sum( apply( clusterTable, 2, max) ) / sum(clusterTable)
    out = results$summary
    out$sampleName = sampleName
    out$depthAdjustment = depthAdjustment
    out$clusterPurity = clusterPurity
    out
  } )
  
  outFinal = rbindlist(pcaTableList)
  outDir = "Dropbox/Documents/statGen/scATACseq/benchmarks/resultsAll/results30"
  outFilename = paste0(sampleName, "DepthNormalization.csv")
  write.csv(outFinal, file.path(outDir, outFilename), row.names = F)
  
}


# plots

plotList = comparePCA(rse, cellTypeName)
outDir = "Dropbox/Documents/statGen/scATACseq/benchmarks/resultsAll/results30"
pdf(file.path(outDir, "DepthNormalization.pdf"), width = 8, height = 3)
g1 = grid.arrange(grobs = plotList, nrow = 1, 
                  left=textGrob(sampleName, rot = 90))
dev.off()
# grid.arrange(grobs = list(g1,g1), nrow = 2) 


plotSinglePCA = function(projection, title, cellType){
  p = qplot(x = projection[ , 1 ], 
            y = projection[ , 2 ], 
            color = cellType)
  p = p + labs(title = title, 
               x="PC1", 
               y="PC2") 
  # p = p + theme(legend.title=element_blank())
  p = p + guides(color=FALSE)
  return(p)
}

comparePCA = function(rse, cellTypeName){
  
  rse = doQC(rse, regionSumCutoff = 5, cellSumCutoffSDs = 3)
  
  cellType = rse[[cellTypeName]]
  X = assay(rse)  
  pca = irlba(t(X), nv = 2)
  projection = t(X) %*% pca$v
  projection_post = projection / colData(rse)$cellSumPostQC
  
  X_pre = t( t(X) / colData(rse)$cellSumPostQC )
  pca_pre = irlba(t(X_pre), nv = 2)
  projection_pre = t(X_pre) %*% pca_pre$v
  
  p1 = plotSinglePCA(projection, "none", cellType)
  
  p2 =  plotSinglePCA(projection_pre, "prePCA", cellType)
  
  p3 =  plotSinglePCA(projection_post, "postPCA", cellType)
  
  return(list(p1, p3, p2))

}

