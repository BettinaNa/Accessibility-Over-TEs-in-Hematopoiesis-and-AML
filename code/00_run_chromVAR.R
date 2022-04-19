#----------------------------------------------------------------
## load dependencies
#----------------------------------------------------------------
suppressMessages(library(chromVAR))
suppressMessages(library(Matrix))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(pheatmap))
suppressMessages(library(data.table))
suppressMessages(library(stringr))

brewer_yes = c('#053061', '#2971B1', '#6AACD0', '#C1DDEB', '#F7F7F7', '#FACDB5', '#E58267', '#BB2933', '#67001F')


run_chromvar=function(peakrds,
                       gsubpattern, 
                       repeatrds, 
                       opname,...){


    # create out_dir
    dir.create(paste0("/results/", opname, "/chromvar/"))
    #----------------------------------------------------------------
    ## Load Countdata
    #----------------------------------------------------------------
    set.seed(2017)

    print("read in count data")
    my_counts_matrix <- readRDS(peakrds)
    rownames(my_counts_matrix) <- paste(my_counts_matrix[,1],my_counts_matrix[,2],my_counts_matrix[,3], sep="_")
    colnames(my_counts_matrix) <- gsub(gsubpattern, "",  colnames(my_counts_matrix))


    ## Remove random chromosomes
    toMatch <- c("random", "alt", "Un", "chrM", "EBV")
    my_counts_matrix <- subset(my_counts_matrix, !(grepl(paste(toMatch, collapse="|"), my_counts_matrix$seqnames)))

    fragment_counts <- makeSummarizedExperimentFromDataFrame(my_counts_matrix)
    assayNames(fragment_counts) <- "counts"

    #----------------------------------------------------------------
    ## add gc content
    #----------------------------------------------------------------
    print("correcting for gc content")
    fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
    counts_filtered <- filterPeaks(fragment_counts,min_fragments_per_peak = 1, non_overlapping = TRUE)
    save(counts_filtered,file=paste0("/results/", opname, "/chromvar/", opname, ".counts_filtered.Rdata"))
    rm(fragment_counts)

    # #----------------------------------------------------------------
    # ## Get motifs and what peaks contain motifs --
    # #----------------------------------------------------------------
    my_annotation_df <- readRDS(repeatrds)
    rownames(my_annotation_df) <- paste(my_annotation_df[,1], my_annotation_df[,2], my_annotation_df[,3], sep="_")
    my_annotation_df <- my_annotation_df[rownames(counts_filtered),]
    rm(my_counts_matrix)
    print("generate annotation database")

    anno_ix <- getAnnotations(as.matrix(my_annotation_df[,4:ncol(my_annotation_df)]), rowRanges = rowRanges(counts_filtered))
    save(anno_ix,file=paste0("/results/", opname, "/chromvar/", opname, ".anno_ix.Rdata"))

    #----------------------------------------------------------------
    ## compute deviation
    #----------------------------------------------------------------
    set.seed(2017)
    print("Computing Deviation")
    dev <- computeDeviations(object = counts_filtered, annotations = anno_ix)
    save(dev,file=paste0("/results/", opname, "/chromvar/", opname, ".devobj.Rdata"))

    z.scores = deviationScores(dev) ## deviation Z-score
    dev.scores = deviations(dev) ## bias corrected deviations

    write.table(z.scores, file=paste0("/results/", opname, "/chromvar/", opname, ".Zscore.txt"), col.names=T, row.names=T, sep="\t", quote=F)
    write.table(dev.scores, file=paste0("/results/", opname, "/chromvar/", opname, ".Deviations.txt"), col.names=T, row.names=T, sep="\t", quote=F)


    #----------------------------------------------------------------
    ## compute variablity
    #----------------------------------------------------------------
    variability <- computeVariability(dev)
    write.table(variability, file=paste0("/results/", opname, "/chromvar/", opname, ".Variability.txt"), col.names=T, row.names=F, sep="\t", quote=F)

}


dirs <- c("hemat", "lsc_hemat", "lsc", "AMLbulk_lsc","AMLbulk") 

for (f in dirs) {
        print(f)
        run_chromvar(peakrds=paste0("/results/",f, "/", f, ".consensus.bed.Binarymat.rds"),
                gsubpattern="",
                repeatrds=paste0("/results/", f,"/",f,".consensus.bed.Binarymat.repeats.rds"),
                opname=f)
}
