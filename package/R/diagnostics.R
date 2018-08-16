# ###Diagnostics-----------------------------------------
# 
# nPeaksVec = seq(2e4,nrow(rse),2e4)
# pValCutoffVec = 2:40
# TSSWeightsList = list(c(0,1), c(1,4), c(1,2), c(1,1.5), c(1,1.25), c(1,1),
#                       c(1.25,1), c(1.5,1), c(2,1), c(4,1), c(1,0))
# DHSWeightsList = list(c(5,1), c(4,1), c(3,1), c(2,1), c(1,1), c(1,2), c(1,3),
#                       c(1,4), c(1,5))
# PCrange = 3:50
# 
# nCores = 24
# cl = makeCluster(nCores)
# clusterEvalQ(cl, library(SummarizedExperiment))
# clusterEvalQ(cl, library(Matrix))
# clusterEvalQ(cl, library(irlba))
# clusterEvalQ(cl, library(data.table))
# clusterExport(cl, list("rse", "TSSWeightsList", "DHSWeightsList", 
#                        "getDestin", "PCrange", "getLogLike", "dmultFast",
#                        "peaksDir", "sampleName", "nPeaksVec", "pValCutoffVec"))
# 
# pValCutoffResults  =  rbindlist(
#   parLapply(cl, pValCutoffVec, function(pValCutoff) {
#     getDestin(rse, pValCutoff = pValCutoff, doLogLikeOG=T)
#   }))
# pValCutoffResults$type = "pValCutoff"
# nPeaksResults  =  rbindlist(
#   parLapply(cl, nPeaksVec, function(myPeak) {
#     getDestin(rse, nPeaks = myPeak, doLogLikeOG=T)
#   }))
# nPeaksResults$type = "nPeaks"
# TSSWeightsResults  =  rbindlist(
#   parLapply(cl, TSSWeightsList, function(TSSWeights) {
#     getDestin(rse, TSSWeights = TSSWeights, doLogLikeOG=T)
#   }))
# TSSWeightsResults$type = "TSSWeights"
# DHSWeightsResults  =  rbindlist(
#   parLapply(cl, DHSWeightsList, function(DHSWeights) {
#     getDestin(rse, DHSWeights = DHSWeights, doLogLikeOG=T)
#   }))
# DHSWeightsResults$type = "DHSWeights"
# nPCResults =  rbindlist(
#   parLapply(cl, PCrange, function(nPCs) {
#     getDestin(rse, PCrange = nPCs, doLogLikeOG=T)
#   }))
# nPCResults$type = "nPCs"
# 
# stopCluster(cl)
# 
# diagnosticsResults = rbind(
#   pValCutoffResults,
#   nPeaksResults,
#   TSSWeightsResults,
#   DHSWeightsResults,
#   nPCResults
# )
# 
# write.csv(diagnosticsResults,
#           file.path(peaksDir, paste0(sampleName, "DestinDiagnostics.csv")),
#           row.names = F)
# 
# optList = list(
#   nPeaks = nPeaksResults[which.max(logLikeOG),]$nPeaksNominal,
#   TSSWeights = c(TSSWeightsResults[which.max(logLikeOG),]$TSSWeight1,
#                  TSSWeightsResults[which.max(logLikeOG),]$TSSWeight2),
#   DHSWeights = c(DHSWeightsResults[which.max(logLikeOG),]$DHSWeight1,
#                  DHSWeightsResults[which.max(logLikeOG),]$DHSWeight2),
#   nPCs = nPCResults[which.max(logLikeOG),]$nPCs
# )
# 
# resultFinalLinear = try(
#   getDestin(rse, PCrange = optList$nPCs, TSSWeight = optList$TSSWeights,
#             DHSWeight = optList$DHSWeights, nPeaks = optList$nPeaks)
# )
# 
# write.csv(resultFinalLinear,
#           file.path(peaksDir, paste0(sampleName, "DestinResultsLinearSummary.csv")),
#           row.names = F)