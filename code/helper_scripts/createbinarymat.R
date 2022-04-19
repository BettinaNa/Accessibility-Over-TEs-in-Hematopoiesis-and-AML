args <- commandArgs(trailingOnly = TRUE)

if (length(args)==3) {
  mergedPeakFile <- as.character(args[1])
  peakFilePattern <- as.character(args[2])
  outputPeakFileName <- as.character(args[3])
  peakFolder <- NULL

} else if (length(args)==4) {
  mergedPeakFile <- as.character(args[1])
  peakFilePattern <- as.character(args[2])
  outputPeakFileName <- as.character(args[3])
  peakFolder <- as.character(args[4])

} else if (length(args)==0) {
  stop("Please provide the following arguments - Merged peak file,  pattern to identify peak files, Binary matrix file name", call.=FALSE)
}

######################################
### load dependencies
######################################
suppressMessages(library(GenomicAlignments))
suppressMessages(library(GenomicRanges))
suppressMessages(library(genomation))
suppressMessages(library(data.table))
suppressMessages(library(stringr))


# create missing variables and run overlap
name <- str_split_fixed(str_split_fixed(mergedPeakFile,"/",4)[,4], "\\.",2)[,1]

outputFolder <-  paste0("/results/", name)
if(!dir.exists(outputFolder)){dir.create(outputFolder)}


if(!is.null(peakFolder)){
  files=dir(path=peakFolder,  pattern = peakFilePattern, full.names=TRUE)
} else if (str_detect(name,"_")){
  # we need to use data from two folders
  peakFolder <- paste0("/data/", str_split_fixed(name, "_",2)[,1] , "_files/")
  files1 = dir(path=peakFolder,  pattern = peakFilePattern, full.names=TRUE)
  peakFolder2 <- paste0("/data/", str_split_fixed(name, "_",2)[,2] , "_files/")
    files2 = dir(path=peakFolder2,  pattern = peakFilePattern, full.names=TRUE)
  files = c(files1, files2)
} else {
    peakFolder <- paste0("/data/", name, "_files/")
    files=dir(path=peakFolder,  pattern = peakFilePattern, full.names=TRUE)
}

query <- suppressMessages(readBed(mergedPeakFile))
queryDF <- data.frame(query)
totalOverlap <- data.frame(seqnames = queryDF$seqnames, start = queryDF$start, end = queryDF$end)
cellName = basename(files)

for (i in 1:length(files)){

  subject = suppressMessages(readBed(files[i]))
  hits = findOverlaps(query, subject)

  hitsDF <- data.frame(hits)
  cellName[i] <- gsub(peakFilePattern, '', cellName[i])
  totalOverlap[hitsDF$queryHits, cellName[i]] <- 1
  totalOverlap[-hitsDF$queryHits, cellName[i]] <- 0

}

outputFile = paste0(outputFolder,"/",outputPeakFileName)
saveRDS(totalOverlap, paste0(outputFile, ".rds"))
