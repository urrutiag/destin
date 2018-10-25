getSpecQuantiles = function(geneAccessibility, nBins, myCellType){
  accessibilityCellType = geneAccessibility[cellType == myCellType,]
  specs = accessibilityCellType$cellTypeSpecificity
  cellQuantiles = quantile(specs[specs!=0], probs = seq(0,1,by=1/nBins))
  accessibilityCellType[, quantile := as.integer(
    cut(cellTypeSpecificity, 
        breaks = cellQuantiles, 
        labels = seq(nBins), 
        include.lowest = T))
    ]
  accessibilityCellType[cellTypeSpecificity == 0, quantile := 0]
  return(accessibilityCellType[])
}

getEWCE = function( geneAccessibility, geneList, fileName ) {
  nHits = sum(geneList %in% geneAccessibility$hgncSymbol)
  empScore= 
    geneAccessibility[hgncSymbol %in% geneList, list(score = sum(cellTypeSpecificity)), by = "cellType"]
  
  # 3) bootstrap score = sum of geneAccessibility[randomGenes,]$cellTypeSpecificity by cell type
  nBoot = 1000
  bootScoreList = rbindlist(lapply(1:nBoot, function(B) {
    bootGeneList = sample(unique(geneAccessibility$hgncSymbol), nHits)
    bootScore= 
      geneAccessibility[hgncSymbol %in% bootGeneList, list(score = sum(cellTypeSpecificity)), by = "cellType"]
    return( bootScore ) 
  }))
  
  #for each cell type
  finalP = rbindlist(lapply(unique(empScore$cellType), function(myCellType) {
    cellScore = empScore[cellType == myCellType]$score
    bootScores = bootScoreList[cellType == myCellType]$score
    
    pValue = mean(bootScores >= cellScore)
    foldChange = cellScore / mean(bootScores)
    sd_from_mean = ( cellScore -mean( bootScores ) ) / sd( bootScores )
    
    myDF = data.frame(cellType = myCellType,
                      pValue = pValue,
                      foldChange = foldChange,
                      sd_from_mean = sd_from_mean,
                      fileName = fileName)
    
    return(myDF)
  }))
  
  #format
  cols = c("foldChange", "sd_from_mean")
  finalP[, (cols) := round(finalP[,mget(cols)],2)]
  
  return(finalP[])
}
