##--------------------------------------------------------------------------------------------
## Objective : For a given query and background file, run LOLA Enrichment
##--------------------------------------------------------------------------------------------

suppressMessages(library(LOLA))
suppressMessages(library(simpleCache))
suppressMessages(library(genomation))

args <- commandArgs(trailingOnly = TRUE)

queryfilePath=as.character(args[1])
unifilePath=as.character(args[2])
lolaregiondbPath=as.character(args[3])
opdirname=as.character(args[4])


## Run LOLA
regionDB = loadRegionDB(lolaregiondbPath)
query = readBed(queryfilePath)
univ = readBed(unifilePath)
res = runLOLA(query, univ, regionDB, cores=1)
writeCombinedEnrichment(res, outFolder=opdirname, includeSplits=TRUE)
