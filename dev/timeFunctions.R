### Functions -----------------------

runAlgorithm = function(pipeResultsDir, sampleName, model, nClusters, 
                        nCells, algorithm, cellData = NULL){
  
  print( paste(sampleName, algorithm) )  
  startTime = Sys.time()
  
  if (algorithm == "Destin"){
    clusterResults = runDestin(pipeResultsDir, model, nClusters, nCells)
  }
  if (algorithm == "Scasat"){
    clusterResults = runScasat(pipeResultsDir, nClusters, nCells)
  }
  if (algorithm == "scABC"){
    clusterResults = runScABC(pipeResultsDir, nClusters, nCells)
  }
  if (algorithm == "chromVar"){
    clusterResults = runChromVar(pipeResultsDir, model, nClusters, nCells)
  }  
  
  runTime = difftime(Sys.time(), startTime, units = "secs")
  print(paste('runTime:', runTime))

  # purity = NA
  # if ( !is.null(cellData) ) {
  #   clusterResults$cellType = 
  #     cellData[ match( clusterResults$cellID, cellData$Run ) ]$source_name
  #   clusterTable = table(clusterResults$cluster, clusterResults$cellType)
  #   purity = sum(apply(clusterTable,1,max))/sum(clusterTable)
  # }

  result = list(sampleName = sampleName,
                algorithm = algorithm,
                runTime = runTime,
                # purity = purity,
                nCells = nrow(clusterResults) )
  

  outDir = "~/Dropbox/Documents/statGen/scATACseq/benchmarks/resultsAll/results31"
  fileName = paste(sampleName, algorithm, 'csv', sep = '.')
  write.csv(data.frame(result),
            file.path( outDir, fileName), 
            row.names = F)
  
  return( result )
}


getPaths = function(pipeResultsDir, nCells){
  
  sampleName = basename(pipeResultsDir)
  bamDir = file.path(pipeResultsDir, "bam")
  peaksDir = file.path(pipeResultsDir, "peaks")
  bedFile = file.path(peaksDir, 
                      paste0(sampleName, "_peaks.blacklist_removed.narrowPeak"))
  bamFiles = dir( bamDir, pattern = "final.bam" )
  bamFiles = bamFiles[!grepl("bai", bamFiles)]
  
  # Preissl exclusions cellType 2
  if (sampleName == "p56"){
    barcodeMap = read.table("/Users/urrutig1/Dropbox/Gene_shared/scATACseq/data/p56_cluster.txt")
    names(barcodeMap) = c("sample", "cellType")
    barcodeMap = data.table(barcodeMap)
    barcodeMap[, cellID := tstrsplit(sample, '_')[[2]] ]

    excludeCells = barcodeMap[ cellType == 2, ]$cellID
    excludeBams = paste0(excludeCells, ".p56.final.bam")
      
    bamFiles = bamFiles[ ! bamFiles %in% excludeBams ]
  }
  
  if ( !is.null(nCells) ) {
    bamFiles = sample(bamFiles, nCells)
  }
  
  paths = list(bamDir = bamDir, 
               bamFiles = bamFiles, 
               bedFile = bedFile)
  
  return(paths)
  
}

getBed = function(bedFile){
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                            qValue = "numeric", peak = "integer")
  bedData = rtracklayer::import(bedFile, format = "BED",
                                extraCols = extraCols_narrowPeak)
  return (bedData)
}


# Destin ---------------------

runDestin = function(pipeResultsDir, model, nClusters, nCells){
  
  print( paste("Running Destin") ) 
  
  paths = getPaths(pipeResultsDir, nCells)
  bedData = getBed(paths$bedFile)
  rse = createRSE(paths$bamDir, paths$bamFiles, bedData)
  
  # sampleName = basename(pipeResultsDir)
  # colData(rse)$cellID = gsub(paste0(".", sampleName, ".final.bam"), "", 
  #                            paths$bamFiles)
    
  rse = annotateRSE(rse, model)
  
  rse = doQC(rse, regionSumCutoff = 5, cellSumCutoffSDs = 3)
  
  sampleName = basename(pipeResultsDir)
  results = destinGrid (rse, nClusters = nClusters, nCores = 7)

  return(results$cluster)
  
}

# scasat ---------------------


runScasat = function(pipeResultsDir, nClusters, nCells){
  
  print( paste("Running scasat") ) 

  paths = getPaths(pipeResultsDir, nCells)
  bedData = getBed(paths$bedFile)
  rse = createRSE(paths$bamDir, paths$bamFiles, bedData)
  
  rse = rse[Matrix::rowSums(assay(rse)) >= 10, ]
  cellSum = Matrix::colSums(assay(rse))
  cutoffs = 2 ^ ( median(log2(cellSum)) + c(-3,3)*mad(log2(cellSum)) )
  rse = rse[ , cellSum > cutoffs[1] & cellSum < cutoffs[2] ]
  
  jacDist = proxy::dist(as(t(assay(rse)), "matrix"), method = "Jaccard")
  kfitScasat = cluster::pam(jacDist, k = nClusters, diss=TRUE)
  
  colData(rse)$cluster = kfitScasat$clustering
  
  return(colData(rse))

}


# scABC ---------------------

runScABC = function(pipeResultsDir, nClusters, nCells){
  
  print( paste("Running scABC") ) 

  paths = getPaths(pipeResultsDir, nCells)

  setwd(paths$bamDir)
  s2 = scABC(paths$bamFiles, paths$bedFile, nClusters = nClusters)
  
  fileNames = names(s2$cluster_assignments)

  clusterResults = data.frame(
    fileName = fileNames,
    cluster = s2$cluster_assignments
  )
  
  return(clusterResults)
  
}


# chromVar --------------

runChromVar = function(pipeResultsDir, model, nClusters, nCells){

  paths = getPaths(pipeResultsDir, nCells)
  peaks = getBed(paths$bedFile)
  bamFiles = paths$bamFiles
  
  bamFiles = bamFiles[!bamFiles %in% 
                        c( "SRR3986179.Corces.final.bam",
                           "SRR3986211.Corces.final.bam")]
  
  colData = data.table::data.table(data.frame(fileName = bamFiles))
  
  setwd(paths$bamDir)
  rse = getCounts(bamFiles, peaks, 
            paired =  TRUE, 
            by_rg = F, 
            format = "bam", 
            colData = colData)
  
  if ( model == "hg19" ) {
    myGenome = BSgenome.Hsapiens.UCSC.hg19
    mySpecies = "Homo sapiens"
  }
  if ( model == "mm10" ) {
    myGenome = BSgenome.Mmusculus.UCSC.mm10
    mySpecies = "Mus musculus"
  }
  
  rse <- addGCBias(rse, genome = myGenome)
  
  rse_filtered <- filterSamples(rse, min_depth = 1500,
                                min_in_peaks = 0.15, shiny = F)
  
  rse_filtered = sort(rse_filtered)
  rse_filtered <- filterPeaks(rse_filtered)
  
  motifs <- getJasparMotifs(species = mySpecies)
  motif_ix <- matchMotifs(motifs, rse_filtered,
                          genome = myGenome)
  
  # computing deviations
  dev <- computeDeviations(object = rse_filtered, 
                           annotations = motif_ix)

  #kmeans
  kfit = kmeans(t(assay(dev)), centers = nClusters, nstart=100)
  
  colData(rse_filtered)$cluster = kfit$cluster
  
  return(colData(rse_filtered))
  
}

# 10x ---------------

cluster10x = function(sampleName, nCells = NULL, nRegions = NULL){
  
  algorithm = "Destin"
  nClusters = 10
  model = "hg19"
  nCores = 7
  # PCrange = 10 # 3:10
  # TSSWeightsList = list(c(1,1)) # list(c(1,2), c(1,1.5), c(1,1.25), c(1,1))
  # DHSWeightsList = list(c(1,1)) #list(c(1,1), c(1,2), c(1,3), c(1,5))
  
  print( paste(sampleName, algorithm) )  
  startTime = Sys.time()
  
  data10xDir = file.path("~/Dropbox/Documents/statGen/scATACseq/10xgenomics",
                         sampleName,
                         "filtered_peak_bc_matrix")
  fileNames = dir(data10xDir, full.names = T)
  
  rse = createRSEfrom10xMatrix(data10xDir)
  
  if ( !is.null(nCells) ) {
    cellSamples = sample(1:ncol(rse), nCells)
    rse = rse[,cellSamples]
  }
  
  if ( !is.null(nRegions) ) {
    regionSamples = sample(1:nrow(rse), nRegions)
    rse = rse[regionSamples,]
  }
  
  rse = annotateRSE(rse, model)
  
  rse = doQC(rse, regionSumCutoff = 5, cellSumCutoffSDs = 3)
  
  clusterResults = destinGrid (rse, nClusters = nClusters, nCores = 7)
  
  runTime = difftime(Sys.time(), startTime, units = "secs")
  print(paste('runTime:', runTime))
  
  result = list(sampleName = sampleName,
                algorithm = algorithm,
                runTime = runTime,
                nCells = nrow(clusterResults$cluster) )
  
  outDir = "~/Dropbox/Documents/statGen/scATACseq/benchmarks/resultsAll/results31"
  fileName = paste(sampleName, algorithm, 'csv', sep = '.')
  write.csv(data.frame(result),
            file.path( outDir, fileName), 
            row.names = F)
}
