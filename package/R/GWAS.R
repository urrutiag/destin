
getEWCE = function( geneAccessibility, geneList, nBoot = 1000 ) {
  nHits = sum(geneList %in% geneAccessibility$hgncSymbol)
  
  if (nHits == 0) { 
    return ("there were no significant snps found in the geneAccessibility data set") 
  }
  
  empScore= 
    geneAccessibility[hgncSymbol %in% geneList, 
                      list(score = sum(cellTypeSpecificity)), by = "cluster"]
  
  # bootstrap score = sum of geneAccessibility[randomGenes,]$cellTypeSpecificity by cluster
  bootScoreList = rbindlist(lapply(1:nBoot, function(B) {
    bootGeneList = sample(unique(geneAccessibility$hgncSymbol), nHits)
    bootScore= 
      geneAccessibility[hgncSymbol %in% bootGeneList, list(score = sum(cellTypeSpecificity)), by = "cluster"]
    return( bootScore ) 
  }))
  
  #for each cluster
  finalP = rbindlist(lapply(unique(empScore$cluster), function(myCluster) {
    cellScore = empScore[cluster == myCluster]$score
    bootScores = bootScoreList[cluster == myCluster]$score
    
    pValue = mean(bootScores >= cellScore)
    foldChange = cellScore / mean(bootScores)
    sd_from_mean = ( cellScore -mean( bootScores) ) / sd( bootScores)
    
    myDF = data.frame(cluster = myCluster,
                      pValue = pValue,
                      foldChange = foldChange,
                      sd_from_mean = sd_from_mean)
    
    return(myDF)
  }))
  
  return(finalP[])
}

annotateMouseToHuman = function(rse){
  
  # MGI symbol (Mouse)
  # mart = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  # conversion = getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), mart=mart)
  load( system.file("mart/mmusculusConversion.Rdata", package = "destin") )
  elementMetadata(rse)$mgi_symbol = 
    conversion$mgi_symbol[
      match(elementMetadata(rse)$feature, 
            conversion$ensembl_gene_id)
      ]
  
  # HGNC (Human) via Mouse-Human homolog
  homData = fread( system.file("mart/HOM_MouseHumanSequence.rpt", package = "destin"), 
                   check.names=TRUE)
  mouseHomData = homData[Common.Organism.Name == "mouse, laboratory", .(Symbol, HomoloGene.ID)]
  elementMetadata(rse)$HomoloGene.ID = 
    mouseHomData$HomoloGene.ID[match(
      elementMetadata(rse)$mgi_symbol,
      mouseHomData$Symbol)]
  humanHomData = homData[Common.Organism.Name == "human", .(Symbol, HomoloGene.ID)]
  elementMetadata(rse)$humanSymbol = 
    humanHomData$Symbol[match(
      elementMetadata(rse)$HomoloGene.ID,
      humanHomData$HomoloGene.ID)]
  # Drop all genes which do not have 1:1 mouse:human orthologs
  rse = rse[!is.na(elementMetadata(rse)$humanSymbol)]
  
  # Entrez ID (Human)
  # mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  # conversion = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene"), mart=mart)
  load( system.file("mart/hsapiensConversion.Rdata", package = "destin") )
  elementMetadata(rse)$entrezID = 
    conversion$entrezgene[match(
      elementMetadata(rse)$humanSymbol,
      conversion$hgnc_symbol)]
  # Drop all genes which do not have entrezID
  rse = rse[!is.na(elementMetadata(rse)$entrezID)]
  
  return( rse )
}



aggregateRSEByGene = function(rse, nCores = NULL){
  
  # create list of aggregated chromatin accessibility by entrez ID 
  # parallel optional
  
  aggregateGene = function(myGene){
    miniRse = rse[rowRanges(rse)$entrezID == myGene]
    outList = list(
      rowData = data.frame(entrezID = myGene,
                           hgncSymbol = rowRanges(miniRse)$humanSymbol[1],
                           mgiSymbol = rowRanges(miniRse)$mgi_symbol[1],
                           nSNPs = nrow(miniRse)),
      assay = t(as(apply(assay(miniRse), 2, sum), "sparseMatrix"))
    )
    return ( outList )
  }
  
  if (!is.null(nCores)){ 
    cl = makeCluster(nCores)
    clusterExport(cl, list("rse", "aggregateGene"), envir = environment() )
    clusterEvalQ(cl, library(SummarizedExperiment, quietly = T))
    seMiniList = parLapply(cl, unique(rowRanges(rse)$entrezID) , 
                           function( myGene ) {
                             aggregateGene(myGene)
                           })
    stopCluster(cl)
  }
  
  if (is.null(nCores)){ 
    seMiniList = lapply(unique(rowRanges(rse)$entrezID) , 
                        function( myGene ) {
                          aggregateGene(myGene)
                        })
  }
  
  # combine list of aggregated chromatin accessibility by gene 
  #https://stackoverflow.com/questions/8843700/creating-sparse-matrix-from-a-list-of-sparse-vectors
  assaysList = lapply(seMiniList, function(myRow) myRow$assay)
  nGenes = length(assaysList)
  nCells = assaysList[[1]]@Dim[2]
  MdfList = lapply(seq_along(assaysList), function(Mindex) {
    i = rep(Mindex, length(assaysList[[Mindex]]@x))
    j = which(diff(assaysList[[Mindex]]@p)==1)
    x = assaysList[[Mindex]]@x
    Mdf = data.frame( i = i, j = j, x = x)
    return(Mdf)
  })  
  Mdf = rbindlist(MdfList)
  assay = sparseMatrix( i = Mdf$i, 
                        j = Mdf$j,
                        x = Mdf$x,
                        dims = c(nGenes, nCells)
  )
  
  #assay = Reduce("rbind", lapply(seMiniList, function(myRow) myRow$assay)), 

  # create summarized experiment
  se = SummarizedExperiment(
    rowData = rbindlist(lapply(seMiniList, function(myRow) myRow$rowData)),
    assay = assay, 
    colData = colData(rse)
  )
  
  return(se)
}


getGeneAccessibility = function(se){
  geneAccessibility = rbindlist(
    lapply(unique(se$cluster), 
           function(myCluster) {
             miniMat = assay(se)[, se$cluster == myCluster]
             peakCount = apply(miniMat,1,sum)/ncol(miniMat)
             outDF = data.frame(hgncSymbol = rowData(se)$hgncSymbol,
                                entrezID = rowData(se)$entrezID,
                                cluster = myCluster,
                                accessibility = peakCount
             )
             return(outDF)
           }
    )
  )
  geneAccessibility[, totalAcrossTypes := sum(accessibility) , by = "entrezID"] 
  geneAccessibility[, cellTypeSpecificity := accessibility/totalAcrossTypes]
  return(geneAccessibility[])
}

getQuantsForMagma = function(geneAccessibility){
  accessbilityClusterQuant = rbindlist(
    lapply( unique( geneAccessibility$cluster ), function(myCluster) {
      getSpecQuantiles(geneAccessibility, nBins=40, myCluster=myCluster)
    }))
  quantsForMagma = dcast(accessbilityClusterQuant, entrezID ~ cluster, value.var = "quantile")
  quantsForMagma
}

getSpecQuantiles = function(geneAccessibility, nBins, myCluster){
  accessbilityCluster = geneAccessibility[cluster == myCluster,]
  specs = accessbilityCluster$cellTypeSpecificity
  cellQuantiles = quantile(specs[specs!=0], probs = seq(0,1,by=1/nBins))
  accessbilityCluster[, quantile := as.integer(
    cut(cellTypeSpecificity, 
        breaks = cellQuantiles, 
        labels = seq(nBins), 
        include.lowest = T))
    ]
  accessbilityCluster[cellTypeSpecificity == 0, quantile := 0]
  return(accessbilityCluster[])
}
