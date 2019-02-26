
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
  countsMat = as(countsMat, "dgCMatrix")
  
  countsMat[countsMat > 1] = 1
  
  colData = data.table::data.table(data.frame(
    fileName = fileNames
  ))
  
  # cellID = gsub(paste0(".",sampleName,".final.bam"),"", fileNames)
  # colData = data.table::data.table(data.frame(
  #   cellID = cellID,
  #   fileName = fileNames
  # ))

  rse = SummarizedExperiment::SummarizedExperiment(assays=list(counts=countsMat),
                             rowRanges=bedData, colData=colData)
  return(rse)
  
}

createRSEfrom10xMatrix = function(data10xDir)
{
  
  barcodes = read.table( file.path(data10xDir, "barcodes.tsv"), sep = "\t")
  
  bedData = rtracklayer::import(file.path(data10xDir, "peaks.bed"), format = "BED")
  
  countsMat = Matrix::readMM(file.path(data10xDir,"matrix.mtx"))
  countsMat = as(countsMat, "dgCMatrix")
  
  countsMat[countsMat > 1] = 1
  
  rse = SummarizedExperiment::SummarizedExperiment(assays=list(counts=countsMat),
                                                   rowRanges=bedData, colData=barcodes)

  # filter out Y chrom
  yIndex = seqnames(rowRanges(rse)) == "chrY"
  rse = rse[!yIndex]
    
  return(rse)
}


annotateRSE = function(rse, model){
  
  if ( model %in% c('hg19', 'hg38') ) {
    filePath = "annotation/encode/DHShg19.Rdata"
  }
  if ( model == 'mm10' ) {
    filePath = "annotation/encode/DHSmm10.Rdata"
  }  

  load( system.file(filePath, package = "destin") )
  
  if  ( model == 'hg38' ) {
    chainPath = system.file("annotation/encode/hg19ToHg38.over.chain", package = "destin")
    chain = import.chain(chainPath)
    DHSrangesLifted = liftOver(DHSranges, chain)
    DHSranges = unlist(DHSrangesLifted)
  }  
  
  regionIndex = findOverlaps(rowRanges(rse), DHSranges, select="first")
  rowRanges(rse)$DHSsum = DHSranges$DHSfrequency[regionIndex]
  
  if ( model == "hg19" ) { 
    data(TSS.human.GRCh37)
    annoFile = TSS.human.GRCh37 
  }
  if ( model == "hg38" ) { 
    data(TSS.human.GRCh38)
    annoFile = TSS.human.GRCh38
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


