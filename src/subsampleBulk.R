#!/usr/bin/env Rscript

#SBATCH --mem=32g
#SBATCH --time=24:00:00

bulkDir="/proj/yuchaojlab/gene/scATACseq/data/CorcesBulk"
subDir="/proj/yuchaojlab/gene/scATACseq/data/CorcesSub"
cellData=file.path(bulkDir, "SraRunTable.txt")

library(data.table)

bamBulkDir = file.path(bulkDir, "bam")
bamSubDir = file.path(subDir, "bam")

bamFilesBulk = dir(bamBulkDir)
bamFilesBulk = bamFilesBulk[!grepl("bai", bamFilesBulk)]

cellTable = fread(cellData)
cellTable = cellTable[Instrument == "NextSeq 500"]
cellIDs = substr(bamFilesBulk, 1, 10)
cellTable = cellTable[Run %in% cellIDs]
cellTypes = unique(cellTable$source_name)


get50subs = function(myCellType) { 

cellIndices = which(cellTable$source_name == myCellType)
sampleIndices = sample(cellIndices, 50, replace = T)

subsampleList = lapply(1:50, function(myIndex) {

  sampleInfo = cellTable[sampleIndices[myIndex]]

  bamFileBulk = bamFilesBulk[grepl(sampleInfo$Run, bamFilesBulk)]
  bamFileSub = paste0(myCellType, myIndex, ".", bamFileBulk)

  subRate = 10^-(2*(runif(1)+ 1))

  subsampleCommand = paste0(
    "samtools view -b -s ", subRate, " ",
    file.path(bamBulkDir, bamFileBulk),  
    " > ", 
    file.path(bamSubDir, bamFileSub)
  )
  system(subsampleCommand)

  indexCommand = paste0(
    "samtools index ",
    file.path(bamSubDir, bamFileSub)
  )
  system(indexCommand)

  wcCommand = paste0(
    "samtools view ",
    file.path(bamSubDir, bamFileSub),
    " | wc -l"
  )
  nReads = system(wcCommand, intern=TRUE)

  sampleInfo$subSampleID = paste0(myCellType, "Sub", myIndex)
  sampleInfo$bamFile = bamFileSub
  sampleInfo$nReads = nReads

  return(sampleInfo)

} )

subSampleTable = rbindlist(subsampleList)

return( subSampleTable ) 

}


subSampleTableAll = rbindlist(
  lapply( cellTypes, function(myCellType) {
    subSampleTable = get50subs(myCellType)
  })
)


write.table(subSampleTableAll, 
    file = file.path(subDir, "SraRunTable.txt"), 
    quote = F,
    row.names = F,
    sep = "\t"
)



