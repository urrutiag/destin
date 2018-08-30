
createRSE = function(bamDir, bamFiles, bedData){
  
  countsList = lapply(bamFiles, function( myBamFile ) {
    bam  = GenomicAlignments::readGAlignments(file.path(bamDir, myBamFile))
    counts = as( GenomicRanges::countOverlaps( bedData, bam ), "dsparseVector")
    out = list( file = myBamFile, counts = counts)
    return(out)
  })
  
  fileNames = sapply(countsList, function(myCounts) {
    myCounts[["file"]]
  })
  
  input = lapply(countsList, function(myCounts) {
    myCounts[["counts"]]
  })
  nRows = length(countsList[[1]]$counts)
  countsMat = Matrix::sparseMatrix(
    x=unlist(lapply(input,slot,"x")),
    i=unlist(lapply(input,slot,"i")),
    p=c(0,cumsum(sapply(input,function(x){length(x@x)}))),
    dims=c(nRows,length(input))
  )
  countsMat[countsMat > 1] = 1
  
  colData = data.table::data.table(data.frame(
    fileName = fileNames
  ))
  
  cellID = gsub(paste0(".",sampleName,".final.bam"),"", fileNames)
  colData = data.table::data.table(data.frame(
    cellID = cellID,
    fileName = fileNames
  ))

  rse = SummarizedExperiment::SummarizedExperiment(assays=list(counts=countsMat),
                             rowRanges=bedData, colData=colData)
  return(rse)
  
}


annotateRSE = function(rse, model){
  
  DHSfile = system.file(
    file.path("annotation/encode", 
              paste0("DHS", model, ".Rdata")), 
    package = "destin")
  
  load(DHSfile)
  regionIndex = findOverlaps(rowRanges(rse), DHSranges, select="first")
  rowRanges(rse)$DHSsum = DHSranges$DHSfrequency[regionIndex]
  
  if ( model == "hg19" ) { 
    data(TSS.human.GRCh37)
    annoFile = TSS.human.GRCh37 
  }
  if ( model == "mm10" ) {
    data(TSS.mouse.GRCm38)
    annoFile = TSS.mouse.GRCm38  
  }
  annotatedPeak = annotatePeakInBatch(rowRanges(rse), 
                                      AnnotationData=annoFile, 
                                      select = "first" )
  rowRanges(rse)$distanceToTSS = annotatedPeak$distancetoFeature
  rowRanges(rse)$feature = annotatedPeak$feature
  rowRanges(rse)$region = "promoter"
  rowRanges(rse)$region[ rowRanges(rse)$distanceToTSS > 1000 | 
                           rowRanges(rse)$distanceToTSS < -2000 ] = 
    "distal element"

  return(rse)                                      
}


