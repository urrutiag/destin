
getDiffAccess = function(rse, clusterMap = NULL, 
                         alternateCellLabel = NULL, nCores = NULL){
  
  # return error if not 
  # all.equal(clusterMap$cellID, rse$cellID)
  
  #INPUT:
  if ( !is.null(clusterMap) ){
    C = clusterMap$cluster
  } else if ( !is.null(alternateCellLabel) ){
    C = alternateCellLabel
  }
  X = assay(rse)
  
  if(!is.null(nCores)){
    cl = makeCluster(nCores)
    clusterEvalQ(cl, library(Matrix))
    clusterExport(cl, list("C", "X"), envir = environment())
  }
  
  pfcList = list()
  for (cellIndex in seq_along(unique(C))) {
    
    Cref = C == unique(C)[cellIndex]
    freqRef = ( Matrix::rowSums(X[,Cref]) + 1) / sum(Cref)
    freqRest = ( Matrix::rowSums(X[,!Cref]) + 1) / sum(!Cref)
    log2FC = log2(freqRef / freqRest)
    
    if(!is.null(nCores)){
      clusterExport(cl, list("Cref"), envir = environment())
      pVals = parSapply(cl, 1:nrow(X), function(regionIndex) {
        mat = table(Cref, X[regionIndex,])
        fish = fisher.test(mat)
        return( fish$p.value )
      } )
    } else {
      pVals = sapply(1:nrow(X), function(regionIndex) {
        mat = table(Cref, X[regionIndex,])
        fish = fisher.test(mat)
        return( fish$p.value )
      } )
    }
    
    pValsCor = p.adjust(pVals, method = "bonferroni")
    out = data.frame(pVals,
                     pValsCor,
                     log2FC)
    names(out) = paste( c("pVals", "pValsCor", "log2FC"),
                        unique(C)[cellIndex], sep = "_")
    pfcList[[cellIndex]] = out
    
  }
  
  if(!is.null(nCores)){
    stopCluster(cl)
  }
  
  pfc = do.call("cbind", pfcList)
  pfc = cbind(pfc)
  colData(rse)$cluster = C
  elementMetadata(rse) = cbind(elementMetadata(rse), pfc)
  
  return(rse)
  
}
