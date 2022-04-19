analysechromvar=function( zscoresdf,dir,qval,annocol1,grouping,opname,dev1,out_path=out_path,...){

    print("Running analysechromvar")
    diff_var <- differentialDeviations3(zscoresdf, grouping, parametric = FALSE)
    colnames(diff_var)[3:4] <- levels(as.factor(grouping))[1:2]
    diff_var <- diff_var[order(diff_var$p_value_adjusted),]
    diff_var$repname <- rownames(diff_var)
    write.table(diff_var, file=paste0(out_path, dir, ".", opname, ".DiffDeviations.with.medians.txt" ), col.names=T, row.names=F, sep="\t", quote=F)
    
    diff_var <- diff_var[,c(1:2,5)]
    write.table(diff_var, file=paste0(out_path, dir, ".",   opname, ".DiffDeviations.txt" ), col.names=T, row.names=F, sep="\t", quote=F)


    top_motifs = subset(diff_var$repname,  diff_var$p_value_adjusted < qval)
    top_devs = as.matrix(zscoresdf[which(rownames(zscoresdf) %in% (top_motifs)), ])
    top_devs[!is.finite(top_devs)] <- NA

    my.breaks <- unique(c(seq(-10, 10, length.out=10)))


    cols <- colorRampPalette(brewer.pal(12, "Paired"))
    mycolors <- cols(length(unique(annocol1$group.V2)))
    names(mycolors) <- unique(annocol1$group.V2)
    mycolors2 <- cols(length(unique(annocol1$group.V3)))
    names(mycolors2) <- unique(annocol1$group.V3)
    mycolors <- list(group.V2=mycolors,group.V3= mycolors2)

    pdf(paste0(out_path, "Heatmap.",dir, ".",   opname, ".DiffDeviations.qval", qval,".pdf" ), onefile=F)
    print(pheatmap(na.omit(top_devs),
    breaks=my.breaks,
    clustering_method="ward.D2",
    clustering_distance_rows="euclidean",
    clustering_distance_cols="euclidean",
    col = brewer_yes,
    annotation_col=annocol1,
    annotation_colors=mycolors,
    cluster_rows=T, cluster_cols=T,
    fontsize_col=5,fontsize_row=3,
    scale="none"))
    dev.off()
}


analysefunction=function( zscoresdf,dir,qval,annocol1,grouping,opname,out_path, fig,...){

    diff_var <- differentialDeviations3(zscoresdf, grouping, parametric = FALSE)
    diff_var <- diff_var[,1:2]
    diff_var <- diff_var[order(diff_var$p_value_adjusted),]
    diff_var$repname <- rownames(diff_var)
    #colnames(diff_var)[3:4] <- levels(as.factor(grouping))
    write.table(diff_var, file=paste0(out_path, dir, ".",   opname, ".DiffDeviations.txt" ), col.names=T, row.names=F, sep="\\t", quote=F)

    top_motifs = subset(diff_var$repname,  diff_var$p_value_adjusted < qval)
    top_devs = as.matrix(zscoresdf[which(rownames(zscoresdf) %in% (top_motifs)), ])
    top_devs[!is.finite(top_devs)] <- NA
    my.breaks <- c(seq(-10, 10, length.out=10))

    cols <- colorRampPalette(brewer.pal(12, "Paired"))
    mycolors <- cols(length(unique(annocol1$group.V2)))
    names(mycolors) <- unique(annocol1$group.V2)
    mycolors2 <- cols(length(unique(annocol1$group.V3)))
    names(mycolors2) <- unique(annocol1$group.V3)
    mycolors3 <- cols(length(unique(annocol1$group.stem)))
    names(mycolors3) <- unique(annocol1$group.stem)
    mycolors5 <- cols(length(unique(annocol1$group.lsc2)))
    names(mycolors5) <- unique(annocol1$group.lsc2)
    mycolors6 <- cols(length(unique(annocol1$group.cyto)))
    names(mycolors6) <- unique(annocol1$group.cyto)
    mycolors7 <- cols(length(unique(annocol1$group.eng)))
    names(mycolors7) <- unique(annocol1$group.eng)
    mycolors8 <- cols(length(unique(annocol1$group.CD34)))
    names(mycolors8) <- unique(annocol1$group.CD34)
    mycolors9 <- cols(length(unique(annocol1$group.CD38)))
    names(mycolors9) <- unique(annocol1$group.CD38)
    mycolors10 <- cols(length(unique(annocol1$group.Gender)))
    names(mycolors10) <- unique(annocol1$group.Gender)
    mycolors11 <- cols(length(unique(annocol1$group.age)))
    names(mycolors11) <- unique(annocol1$group.age)
    mycolors12 <- cols(length(unique(annocol1$group.BulkName)))
    names(mycolors12) <- unique(annocol1$group.BulkName)
    mycolors13 <- cols(length(unique(annocol1$group.CD34CD38)))
    names(mycolors13) <- unique(annocol1$group.CD34CD38)

    mycolors <- list(group.V2=mycolors,group.V3= mycolors2,group.stem=mycolors3,
                     group.lsc2= mycolors5,group.cyto=mycolors6,group.eng=mycolors7,
                     group.CD34=mycolors8, group.CD38=mycolors9, group.Gender=mycolors10,
                     group.age=mycolors11, group.BulkName=mycolors12, group.CD34CD38=mycolors13)

    pdf(paste0(out_path, "Fig",fig,"_Heatmap.",dir, ".",   opname, ".DiffDeviations.qval", qval,".pdf" ), onefile=F)
    pheatmap(na.omit(top_devs),
    breaks=my.breaks,
    clustering_method="ward.D2",
    clustering_distance_rows="euclidean",
    clustering_distance_cols="euclidean",
    col = brewer_yes,
    annotation_col=annocol1,
    annotation_colors=mycolors,
    cluster_rows=T, cluster_cols=T,
    fontsize_col=5,fontsize_row=3,
    scale="none")
    dev.off()
}

# compute LSCTE121 Z-score
LSCTE121_Zscore_calculation <- function(opname, path_counts_filtered, repeat_consensus_bed_path, out_dir){
    #load filtered counts
    load(paste0(path_counts_filtered, opname, ".counts_filtered.Rdata"))

   # load repeat annotations
    print("load annotations")
    repeat_annotation_files <- c(paste0(repeat_consensus_bed_path,"LSC_enriched.bed.sorted"), paste0(repeat_consensus_bed_path,"LSC_depleted.bed.sorted")) 
    anno_rep <- getAnnotations(repeat_annotation_files,rowRanges = rowRanges(counts_filtered))

    #set seed 
    set.seed(30)

    #compute deviation
    print("compute deviation")
    dev_rep <- computeDeviations(object = counts_filtered, annotations = anno_rep)

    save(dev_rep,file=paste0(out_dir, opname, "_rep_deviations.RData"))

    # calculate LSCTE121 Z-score
    newcoldata<-colData(dev_rep)
    zscores_rep<-assays(dev_rep)$z
    newcoldata$LSC_enriched_zscores <- zscores_rep[1,]
    newcoldata$LSC_depleted_zscores <- zscores_rep[2,]
    newcoldata$LSC_common_zscore <- (0.5*newcoldata$LSC_enriched_zscores - 0.5*newcoldata$LSC_depleted_zscores)/sqrt((0.5^2+0.5^2))

    # save file
    write.csv(newcoldata, file=paste0(out_dir,opname,"_LSCTE121_Zscores.csv"))

    return(newcoldata)
}


analysechromvar_medians=function( zscoresdf,dir,grouping,opname,...){
    diff_var <- differentialDeviations5(zscoresdf, grouping, parametric = FALSE)
    diff_var$repname <- rownames(diff_var)
    write.table(diff_var, file=paste0(out_path, dir, ".",   opname, ".DiffDeviations.with.medians.txt" ), col.names=T, row.names=F, sep="\t", quote=F)
}

