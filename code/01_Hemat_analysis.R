#----------------------------------------------------------------
## Hematopoietic cell type analysis
#----------------------------------------------------------------
## load dependencies
#----------------------------------------------------------------
rm(list = ls())
suppressMessages(library(chromVAR))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))
suppressMessages(library(UpSetR))

source("/code/helper_scripts/chromvar.helperfunctions.R")
source("/code/helper_scripts/chromvar_analyse_functions.R")

'%!in%' <- function(a,b) ! a %in% b

brewer_yes = c('#053061', '#2971B1', '#6AACD0', '#C1DDEB', '#F7F7F7', '#FACDB5', '#E58267', '#BB2933', '#67001F')

# define out_path
out_path <- "/results/manuscript_figures/Hemat/"
if(!dir.exists(out_path)){dir.create(out_path, recursive = TRUE)}

#----------------------------------------------------------------
## Read in metadata
#----------------------------------------------------------------
# this metadata has the ST-HSCs marked as stem populations
metadata <- read.table("/data/metadata/metadata_ltst_stem.txt", header=F,sep="\t", stringsAsFactors=F)
metadata <- unique(metadata)
rownames(metadata) <- str_split_fixed(metadata[,1],"\\.",2)[,1]

metadata$Stem <- ifelse(metadata$V3 == "Hemat_diff","Mature","Primitive")

#------------------------------------------------------------------------------------
## hemat
#------------------------------------------------------------------------------------
zscore_hemat <- read.table("/results/hemat/chromvar/hemat.Zscore.txt",header=T, check.names=F,sep="\t", stringsAsFactors=F)
annocol_hemat <- data.frame(stringsAsFactors=F, group=metadata[colnames(zscore_hemat) ,2:3])
load("/results/hemat/chromvar/hemat.devobj.Rdata")

#---------------------------------------------
## Stem and Progenitor VS Differentiated 
#---------------------------------------------
annocol_hemat2 <- annocol_hemat
annocol_hemat2$group.V3 <- ifelse(annocol_hemat2$group.V3 %in% c("Hemat_prog","Hemat_stem"),"Primitive","Mature" )

analysechromvar(zscoresdf=zscore_hemat,
                dir="Hemat_TEonly",qval=0.01,
                annocol1=annocol_hemat2,
                grouping=as.factor(annocol_hemat2$group.V3),
                opname="Hemat_stemprog_vs_diff",dev1=dev,
                out_path=out_path)
rm(dev)


#-------------------------------------------------------
## Re-plot heatmap for Hemat_stemprog_vs_diff
## Fig 1a
#-------------------------------------------------------
# adjust colors to correct coloring: "#F36F60", "#F79C76", "#A7D1D5"

cols <- colorRampPalette(brewer.pal(12, "Paired"))
mycolors <- cols(length(unique(annocol_hemat$group.V2)))
names(mycolors) <- unique(annocol_hemat$group.V2)
mycolors2 <- c("#F36F60", "#F79C76", "#A7D1D5")
names(mycolors2) <- unique(annocol_hemat$group.V3)
mycolors <- list(group.V2=mycolors,group.V3= mycolors2)

my.breaks <- c(seq(-10, 10, length.out=10))

pdf(paste0(out_path,"/Fig1a_Heatmap.Hemat_TEonly.Hemat_stemprog_vs_diff.allRepeats.pdf"), onefile=F)
pheatmap(na.omit(zscore_hemat),
breaks=my.breaks,
clustering_method="ward.D2",
clustering_distance_rows="euclidean",
clustering_distance_cols="euclidean",
col = brewer_yes,
annotation_col=annocol_hemat,
annotation_colors=mycolors,
cluster_rows=T, cluster_cols=T,
fontsize_col=5,fontsize_row=3,
show_rownames=F,show_colnames=F,
scale="none")
dev.off()


#-------------------------------------------------------
## Plot volcano plot for Hemat_stemprog_vs_diff
#-------------------------------------------------------
# Suppl Fig 1d
#----------------------------
# adjusted to q-value: 0.01
Hemat_stemprog_vs_diff.diffdeviations <- read.table(paste0(out_path,"/Hemat_TEonly.Hemat_stemprog_vs_diff.DiffDeviations.with.medians.txt"), header=T, sep="\t", stringsAsFactors=F)
Hemat_stemprog_vs_diff.diffdeviations$mediandiff <- Hemat_stemprog_vs_diff.diffdeviations$Primitive-Hemat_stemprog_vs_diff.diffdeviations$Mature
Hemat_stemprog_vs_diff.diffdeviations$neglog10qval <- -log10(Hemat_stemprog_vs_diff.diffdeviations$p_value_adjusted)
Hemat_stemprog_vs_diff.diffdeviations$type <- ifelse(Hemat_stemprog_vs_diff.diffdeviations$p_value_adjusted < 0.01, "sig", "nonsig")
Hemat_stemprog_vs_diff.diffdeviations$type2 <- ifelse(Hemat_stemprog_vs_diff.diffdeviations$type=="sig" & Hemat_stemprog_vs_diff.diffdeviations$mediandiff>0, "Primitive",
                                                 ifelse(Hemat_stemprog_vs_diff.diffdeviations$type=="sig" & Hemat_stemprog_vs_diff.diffdeviations$mediandiff<0,"Mature", "nonsig"))

pdf(paste0(out_path,"/SupplFig1d_Hemat_stemprog_vs_diff.volcanoplot.pdf"), useDingbats=F)
ggplot(data=Hemat_stemprog_vs_diff.diffdeviations, aes(x=mediandiff, y=neglog10qval)) +
  geom_point(aes(color=type2)) +  
  theme(text = element_text(size=12)) + 
  labs(x = "Median difference (Primitive - Mature)" , y="Negative log q-value") + 
  scale_colour_manual(values=c(Primitive="#CC0033",nonsig="#000000", Mature="#1e71aa")) +
  geom_hline(yintercept=2,linetype="dashed", size=0.5, col="grey") +
  xlim(-15, 10) +   theme_classic()
dev.off()

#---------------------------
# Fig 1b
# family overview
#---------------------------
# load TE family annotation
repeat_metadata <- read.table("/data/metadata/repeat_metadata.txt", header=T, sep="\t", stringsAsFactors=F)
repfam <- read.table("/data/metadata/repeat_fam_mapping", header=F , sep="\t", stringsAsFactors=F)
repeat_metadata <- merge(repeat_metadata, repfam[,1:2], by.x="subfam",by.y="V2", all.x=T)
repeat_metadata[is.na(repeat_metadata)] <- "UNKNOWN"

# merge annotation with TE family annotation
repfam_sig_hemat <- merge(Hemat_stemprog_vs_diff.diffdeviations, repeat_metadata, by="repname", all.x=TRUE, all.y=FALSE)

repfam_sig_hemat$hemat_sig <- ifelse(repfam_sig_hemat$type == "sig", 1,0)
repfam_sig_hemat$hemat_sigstem <- ifelse(repfam_sig_hemat$type2 == "Primitive", 1,0)
repfam_sig_hemat$hemat_signonstem <- ifelse(repfam_sig_hemat$type2 == "Mature", 1,0)

# make summary table to determine which families do not show significant subfamilies
summary_families_hemat <- plyr::ddply(repfam_sig_hemat, "V1", summarise, count=length(type), sig=sum(hemat_sig), primitive =sum(hemat_sigstem), mature=sum(hemat_signonstem), .drop=FALSE)
write.csv(summary_families_hemat, paste0(out_path,"/Hemat_repfam_overview.csv"))

# don't plot families without significant subfamilies
not <- summary_families_hemat[summary_families_hemat$sig == 0,][,1]

repfam_sig_hemat <- repfam_sig_hemat[!(repfam_sig_hemat$V1 %in% not),]
repfam_sig_hemat$famSort <- factor(repfam_sig_hemat$V1, )


pdf(paste0(out_path,"/Fig1b_Hemat_TEonly_repfam_overview_jitter_type2_onlysig_q0.01_TEonly.pdf" ), onefile=F, useDingbats=FALSE)
ggplot(repfam_sig_hemat,aes(x=famSort,y=mediandiff))+
    geom_jitter(aes(colour=type2), width=0.1)+
    coord_flip()+
    theme_bw() + 
    scale_x_discrete(limits = rev(levels(repfam_sig_hemat$famSort))) + 
    scale_colour_manual(values=c(Primitive="#CC0033",nonsig="#000000", Mature="#1e71aa"))+ 
    labs(x = "TE family" , y="Median difference (Primitive - Mature)") 
dev.off()


##------------------------------------------------------
# Find Ery lineage specific TEs
##------------------------------------------------------
# Suppl Fig 1b
#-----------------------
annocol_hemat3 <- annocol_hemat
annocol_hemat3$group.V4 <- ifelse(annocol_hemat3$group.V2 %in% c("MEP","Erythroid Precustor"),
                                          "Ery","non_Ery" )

analysechromvar(zscoresdf=zscore_hemat,
                dir="SupplFig1b",qval=0.001,
                annocol1=annocol_hemat3,
                grouping=as.factor(annocol_hemat3$group.V4),
                opname="Hemat_EryLin_nonEry",
                dev1=dev, 
                out_path=out_path)


Hemat_EryLin_nonEry.diffdeviations <- read.table(paste0(out_path,"/SupplFig1b.Hemat_EryLin_nonEry.DiffDeviations.with.medians.txt"), header=T, sep="\t", stringsAsFactors=F)
Hemat_EryLin_nonEry.diffdeviations.reps <- subset(Hemat_EryLin_nonEry.diffdeviations, Hemat_EryLin_nonEry.diffdeviations$p_value_adjusted <0.001)

#--------------------------
# Suppl Fig 1c
#--------------------------
zscore_hemat_ery <- subset(zscore_hemat,rownames(zscore_hemat) %!in% Hemat_EryLin_nonEry.diffdeviations.reps$repname)
cols <- colorRampPalette(brewer.pal(12, "Paired"))

mycolors <- cols(length(unique(annocol_hemat3$group.V2)))
names(mycolors) <- unique(annocol_hemat3$group.V2)
mycolors2 <- c("#F36F60", "#F79C76", "#A7D1D5")
names(mycolors2) <- unique(annocol_hemat3$group.V3)
mycolors3 <- cols(length(unique(annocol_hemat3$group.V4)))
names(mycolors3) <- unique(annocol_hemat3$group.V4)
mycolors <- list(group.V2=mycolors,group.V3= mycolors2, group.V4=mycolors3)

my.breaks <- c(seq(-10, 10, length.out=10))


pdf(paste0(out_path,"/SupplFig1c_Heatmap.Hemat_TEonly.Hemat_stem_prog_diff.allRepeats_woEryLin28.pdf"), onefile=F)
pheatmap(na.omit(zscore_hemat_ery),
breaks=my.breaks,
clustering_method="ward.D2",
clustering_distance_rows="euclidean",
clustering_distance_cols="euclidean",
col = brewer_yes,
annotation_col=annocol_hemat3,
annotation_colors=mycolors,
cluster_rows=T, cluster_cols=T,
fontsize_col=5,fontsize_row=3,
show_rownames=F,show_colnames=F,
scale="none")
dev.off()



#---------------------------------------------
# Plot boxplots for top TEs
#---------------------------------------------
##---------------------
# Fig 1c
##---------------------
zscore_anno <- t(zscore_hemat)
colnames(zscore_anno) <- str_split_fixed(colnames(zscore_anno),"_",2)[,2]
zscore_anno <- merge(annocol_hemat, zscore_anno, by=0, all=TRUE)
rownames(zscore_anno) <- zscore_anno[,1]
zscore_anno <- zscore_anno[,-1]
zscore_anno <- zscore_anno[order(zscore_anno$group.V3),]

zscore_anno$cell_type <-factor(zscore_anno$group.V2,
    levels = c("LT-HSC","ST-HSC", "MLP", "CMP", "GMP", "MEP", "T Cell", "B Cell","Erythroid Precustor", "Monocyte", "Granulocyte", "Dendritic Cell", "NK Cell"),ordered = TRUE)
zscore_anno$grouping <- ifelse(zscore_anno$group.V3 == "Hemat_stem", "Stem", 
                          ifelse(zscore_anno$group.V3 == "Hemat_prog", "Prog", "Mature"))
zscore_anno$grouping <- factor(zscore_anno$grouping, levels=c("Stem","Prog","Mature"), ordered=TRUE)

# top TEs
non_stem <- c("AluJb","AluSx","MIR3","MIR")
stem <- c("LTR33","LTR16E1","LTR39","MER61A")

plotlist <- list()
for (i in c(stem, non_stem)){
    # use TE family annotation in plot name
    fam <- subset(repeat_metadata, repeat_metadata$subfam == eval(i))[,"V1"]
    plotlist[[eval(i)]] <- ggplot(zscore_anno, aes_string(x="grouping", y=eval(i), fill="grouping")) + 
                            geom_boxplot() + 
                            geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill="black") + 
                            theme_bw() + 
                            scale_fill_manual(values=c("#F36F60", "#F79C76", "#A7D1D5"))+
                            labs(subtitle=paste0(eval(i)," - ",fam)) + 
                            xlab("Group") +
                            ylab("Z-score")+
                            theme(legend.position = "none")
}

pdf(paste0(out_path, "Fig1c_Stem_prog_diff_repfam_overview_with_dotplot_all_final.pdf" ), onefile=F, useDingbats=FALSE)
plot_grid(plotlist = plotlist, nrow=2, ncol=4)
dev.off()

#----------------------------------
# Supplementary figure 1e
# plot previous plot by cell type
#----------------------------------
plotlist <- list()
for (i in c(stem, non_stem)){
    # use TE family annotation in plot name
    fam <- subset(repeat_metadata, repeat_metadata$subfam == eval(i))[,"V1"]
    plotlist[[eval(i)]] <- ggplot(zscore_anno, aes_string(x="cell_type", y=eval(i), fill="grouping")) + 
                            geom_boxplot() + 
                            geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill="black") + 
                            theme_bw() + 
                            coord_flip() +
                            scale_fill_manual(values=c("#F36F60", "#F79C76", "#A7D1D5"))+
                            labs(subtitle=paste0(eval(i)," - ",fam)) + 
                            xlab("Group") +
                            ylab("Z-score")+
                            theme(legend.position = "none")
}

pdf(paste0(out_path, "SupplFig1e_Stem_prog_diff_repfam_overview_with_dotplot_by_cell_type.pdf" ), onefile=F, useDingbats=FALSE, height=8, width=14)
plot_grid(plotlist = plotlist, nrow=2, ncol=4)
dev.off()




#----------------------------------
# Supplementary figure 1f
#----------------------------------

# peak number for supplementary table
peak_numbers <- read.delim("/data/metadata/peak_numbers_hemat_bed.txt", stringsAsFactors=F, header=F, sep=" ")
peak_numbers$names <- gsub(".narrowPeak","", peak_numbers$V5)

metadata$names <- gsub("_peaks", "", gsub("\\..*", "", rownames(metadata)))
# merge with metadata
meta <- merge(peak_numbers[,c("V4","V5","names")], metadata, by="names", all.x=T,all.y=F)
# plot split by cell types
meta$V2 <- factor(meta$V2, levels=c("LT-HSC","ST-HSC", "MLP", "CMP", "GMP", "MEP", "T Cell", "B Cell","Erythroid Precustor", "Monocyte", "Granulocyte", "Dendritic Cell", "NK Cell"))
meta$V3 <- factor(meta$V3, levels=c("Hemat_stem", "Hemat_prog", "Hemat_diff"))

pdf(paste0(out_path, "SupplFig1f_Peak_number_by_cell_types_with_dotplot.pdf" ), onefile=F, useDingbats=FALSE, heigh=7, width=5)
ggplot(meta, aes(x=V2, y=V4)) + 
        geom_boxplot(aes(fill=meta$V3)) + 
        geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill="black") + 
        theme_bw() + 
        scale_fill_manual(values=c("#F36F60", "#F79C76", "#A7D1D5"))+
        xlab("Cell type") +
        ylab("# of peaks")+
        theme(legend.position = "none", axis.text.x=element_text(angle=90, hjust=1))
dev.off()



#----------------------------------
# Suppl Fig 2d
#----------------------------------

subfamilies <- c("MamGypLTR1b", "MamGypsy2-LTR", "MER57B1", "LTR67B", "MER57B2")
colnames(zscore_anno) <- gsub("-","_", colnames(zscore_anno))
zscore_anno1 <- subset(zscore_anno, zscore_anno$cell_type != "ST-HSC")
# loop through plotting
plotlist2 <- list()
for (i in subfamilies){
    # use TE family annotation in plot name
    fam <- subset(repeat_metadata, repeat_metadata$subfam == eval(i))[,"V1"]
    i <- gsub("-","_", i)
    plotlist2[[eval(i)]] <- ggplot(zscore_anno1, aes_string(x="grouping", y=eval(i), fill="grouping")) + 
                            geom_boxplot() + 
                            geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill="black") + 
                            theme_bw() + 
                            scale_fill_manual(values=c("#F36F60", "#F79C76", "#A7D1D5"))+
                            labs(subtitle=paste0(eval(i)," - ",fam)) + 
                            xlab("Group") +
                            ylab("Z-score")+
                            theme(legend.position = "none")
}

pdf(paste0(out_path, "SupplFig2d_LT_HSC_prog_diff_repfam_SupplFig2d.pdf" ), onefile=F, useDingbats=FALSE, height=3, width=9)
plot_grid(plotlist = plotlist2, nrow=1, ncol=5)
dev.off()

 

#----------------------------------
# Suppl Fig 2a
#----------------------------------
#---------------------------------------------
## LT-HSCs only VS mature 
#---------------------------------------------
annocol_hemat2 <- annocol_hemat
annocol_hemat2$group.V3 <- ifelse(annocol_hemat2$group.V3 %in% c("Hemat_prog","Hemat_stem"),"Primitive","Mature" )
annocol_hemat2 <- subset(annocol_hemat2,annocol_hemat2$group.V2 == "LT-HSC" |  annocol_hemat2$group.V3 == "Mature")


zscore_hemat1 <- zscore_hemat[, rownames(annocol_hemat2)]
analysechromvar(zscoresdf=zscore_hemat1,
                dir="Hemat_TEonly",qval=0.01,
                annocol1=annocol_hemat2,
                grouping=as.factor(annocol_hemat2$group.V3),
                opname="LT.HSC_vs_mature",dev1=dev,
                out_path=out_path)

# subset to <0.01
LT.HSC_vs_mature.diffdeviations <- read.table(paste0(out_path,"/Hemat_TEonly.LT.HSC_vs_mature.DiffDeviations.with.medians.txt"), header=T, sep="\t", stringsAsFactors=F)
lt_mature <- subset(LT.HSC_vs_mature.diffdeviations, LT.HSC_vs_mature.diffdeviations$p_value_adjusted <0.01 )

annocol_hemat1 <- subset(annocol_hemat, annocol_hemat$group.V3 %!in% "Hemat_prog" & annocol_hemat$group.V2 %!in% "ST-HSC")

cols <- colorRampPalette(brewer.pal(12, "Paired"))
mycolors <- cols(length(unique(annocol_hemat$group.V2)))
names(mycolors) <- unique(annocol_hemat$group.V2)
mycolors <- mycolors[which(names(mycolors) %in% unique(annocol_hemat1$group.V2))]

mycolors2 <- c("#F36F60", "#F79C76", "#A7D1D5")
names(mycolors2) <- unique(annocol_hemat$group.V3)
mycolors2 <-  mycolors2[which(names(mycolors2) %in% unique(annocol_hemat1$group.V3))]

mycolors <- list(group.V2=mycolors,group.V3= mycolors2)

df <- na.omit(zscore_hemat[lt_mature$repname,rownames(annocol_hemat1)])
my.breaks <- c(seq(-10, 10, length.out=10))

pdf(paste0(out_path,"/SupplFig2a_Heatmap.LT-HSCs.vs.Mature.populations.pdf"), onefile=F)
pheatmap(df,
breaks=my.breaks,
clustering_method="ward.D2",
clustering_distance_rows="euclidean",
clustering_distance_cols="euclidean",
col = brewer_yes,
annotation_col=annocol_hemat1,
annotation_colors=mycolors,
cluster_rows=T, cluster_cols=T,
fontsize_col=5,fontsize_row=3,
show_rownames=F,show_colnames=F,
scale="none")
dev.off()


#----------------------------------
# Suppl Fig 2b
#----------------------------------
#---------------------------------------------
## Progenitors VS mature 
#---------------------------------------------

annocol_hemat2 <- annocol_hemat
annocol_hemat2 <- subset(annocol_hemat2,annocol_hemat2$group.V3 != "Hemat_stem")
annocol_hemat2$group.V3 <- ifelse(annocol_hemat2$group.V3 %in% c("Hemat_prog","Hemat_stem"),"Primitive","Mature" )


zscore_hemat1 <- zscore_hemat[, rownames(annocol_hemat2)]
analysechromvar(zscoresdf=zscore_hemat1,
                dir="Hemat_TEonly",qval=0.01,
                annocol1=annocol_hemat2,
                grouping=as.factor(annocol_hemat2$group.V3),
                opname="Prog_vs_mature",dev1=dev,
                out_path=out_path)

# subset to <0.01
Prog_vs_mature.diffdeviations <- read.table(paste0(out_path,"/Hemat_TEonly.Prog_vs_mature.DiffDeviations.with.medians.txt"), header=T, sep="\t", stringsAsFactors=F)
prog_mature <- subset(Prog_vs_mature.diffdeviations, Prog_vs_mature.diffdeviations$p_value_adjusted <0.01 )


annocol_hemat1 <- subset(annocol_hemat, annocol_hemat$group.V3 %!in% "Hemat_stem" )

cols <- colorRampPalette(brewer.pal(12, "Paired"))
mycolors <- cols(length(unique(annocol_hemat$group.V2)))
names(mycolors) <- unique(annocol_hemat$group.V2)
mycolors <- mycolors[which(names(mycolors) %in% unique(annocol_hemat1$group.V2))]

mycolors2 <- c("#F36F60", "#F79C76", "#A7D1D5")
names(mycolors2) <- unique(annocol_hemat$group.V3)
mycolors2 <-  mycolors2[which(names(mycolors2) %in% unique(annocol_hemat1$group.V3))]

mycolors <- list(group.V2=mycolors,group.V3= mycolors2)

df <- na.omit(zscore_hemat[prog_mature$repname,rownames(annocol_hemat1)])
my.breaks <- c(seq(-10, 10, length.out=10))


pdf(paste0(out_path,"/SupplFig2b_Heatmap.Progenitors.noST_vs_Mature.pdf"), onefile=F)
pheatmap(df,
breaks=my.breaks,
clustering_method="ward.D2",
clustering_distance_rows="euclidean",
clustering_distance_cols="euclidean",
col = brewer_yes,
annotation_col=annocol_hemat1,
annotation_colors=mycolors,
cluster_rows=T, cluster_cols=T,
fontsize_col=5,fontsize_row=3,
show_rownames=F,show_colnames=F,
scale="none")
dev.off()


##---------------------------------
## Supplementary Figure 2c
##---------------------------------
# upset plot - LT-HSC vs Mature, Progenitors vs Mature and Primitive TEs
Hemat_stemprog_vs_diff.diffdeviations.Primitive_reps <- subset(Hemat_stemprog_vs_diff.diffdeviations, Hemat_stemprog_vs_diff.diffdeviations$type2=="Primitive")
prog_mature.Prog <- subset(Prog_vs_mature.diffdeviations, Prog_vs_mature.diffdeviations$p_value_adjusted <0.01 & Prog_vs_mature.diffdeviations$Primitive > 0)
lt_mature.lt <- subset(LT.HSC_vs_mature.diffdeviations, LT.HSC_vs_mature.diffdeviations$p_value_adjusted <0.01 & LT.HSC_vs_mature.diffdeviations$Primitive > 0)


pdf(paste0(out_path, "SupplFig2c_UpsetR_plot_LT-HSC_Prog_Primitive_TEs.pdf"), onefile=F)
upset(fromList(list(Primitive=Hemat_stemprog_vs_diff.diffdeviations.Primitive_reps$repname, `LT-HSC`=lt_mature.lt$repname, Progenitors=prog_mature.Prog$repname)),order.by = "freq") 
dev.off()






