#-------------------------------------------------------
## AML bulk analysis
#-------------------------------------------------------
rm(list = ls())
suppressMessages(library(chromVAR))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tibble))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))

source("/code/helper_scripts/chromvar.helperfunctions.R")
source("/code/helper_scripts/chromvar_analyse_functions.R")

out_path <- "/results/manuscript_figures/AMLbulk/"
if(!dir.exists(out_path)){dir.create(out_path, recursive = TRUE)}

brewer_yes = c('#053061', '#2971B1', '#6AACD0', '#C1DDEB', '#F7F7F7', '#FACDB5', '#E58267', '#BB2933', '#67001F')

#----------------------------------------------------------------
## Read in metadata
#----------------------------------------------------------------
metadata <- read.csv("/data/metadata/curated_metadata_lsc_hemat_new.csv", row.names=2, stringsAsFactors=F)
metadata <- metadata[,-c(1,2)]
metadata$CD34CD38 <- ifelse(metadata$CD34==1 & metadata$CD38==0, "CD34pos.CD38neg",
                        ifelse(metadata$CD34==1 & metadata$CD38==1, "CD34pos.CD38pos",
                        ifelse(metadata$CD34==0 & metadata$CD38==1, "CD34neg.CD38pos",
                        ifelse(metadata$CD34==0 & metadata$CD38==0, "CD34neg.CD38neg",NA))))

VennTable <- read.table("/results/manuscript_figures/LSC/VennTable_TEonly.txt", stringsAsFactors = F)
LSC <- subset(VennTable, VennTable$group %in% c("LSC+", "LSC-"))[,1]


#----------------------------------------------------------------
## LSCs and bulk AML samples
#----------------------------------------------------------------
zscore_AMLbulk_lsc <- read.table("/results/AMLbulk_lsc/chromvar/AMLbulk_lsc.Zscore.txt",header=T, sep="\t", stringsAsFactors=F)
colnames(zscore_AMLbulk_lsc) <- gsub("X", "", colnames(zscore_AMLbulk_lsc))

# calculate LSCTE121 Z-score
LSCTE121_Zscore <- LSCTE121_Zscore_calculation(opname="AMLbulk_lsc", path_counts_filtered="/results/AMLbulk_lsc/chromvar/", repeat_consensus_bed_path="/data/metadata/", out_dir=out_path)

# make heatmap annotation
# annocol <- merge(metadata[colnames(zscore_AMLbulk_lsc),c("LSC17_score", "LSC_high_low","V2")], as.data.frame(LSCTE121_Zscore[,"LSC_common_zscore"]), by=0) %>%
#             column_to_rownames("Row.names")
annocol <- merge(metadata[colnames(zscore_AMLbulk_lsc),c("LSC_high_low","V2")], as.data.frame(LSCTE121_Zscore[,"LSC_common_zscore"]), by=0) %>%
            column_to_rownames("Row.names")

colnames(annocol)[3] <- "LSCTE121_Z_score"
annocol$V2 <- ifelse(annocol$V2=="LSC.Bulk", "AMLbulk", annocol$V2)

# plot heatmap for LSC repeats
my.breaks <- c(seq(-10, 10, length.out=10))

cols <- colorRampPalette(brewer.pal(12, "Paired"))

mycolors <- cols(length(unique(annocol$V2)))
names(mycolors) <- unique(annocol$V2)
mycolors1 <- cols(length(unique(annocol$LSC_high_low)))
names(mycolors1) <- unique(annocol$LSC_high_low)

mycolors <- list(V2=mycolors,LSC_high_low=mycolors1)

##----------------------
# Figure 4a
##----------------------
# plot for LSC repeats

pdf(paste0(out_path, "Fig4a_Heatmap.AMLbulk_LSCs.LSCrepeats_q0.01.pdf" ), onefile=F)
pheatmap(zscore_AMLbulk_lsc[LSC,],
breaks=my.breaks,
clustering_method="ward.D2",
clustering_distance_rows="euclidean",
clustering_distance_cols="euclidean",
col = brewer_yes,
annotation_col=annocol,
annotation_colors=mycolors,
cluster_rows=T, cluster_cols=T,
fontsize_col=5,fontsize_row=3,
fontsize =5,
show_rownames=F,
show_colnames=F,
scale="none")
dev.off()

#----------------------------------------------------------------
## AMLbulk only
#----------------------------------------------------------------
zscore_AMLbulk <- read.table("/results/AMLbulk/chromvar/AMLbulk.Zscore.txt",header=T, sep="\t", stringsAsFactors=F)
colnames(zscore_AMLbulk) <- gsub("X", "", colnames(zscore_AMLbulk))

# calculate LSCTE121 Z-score
LSCTE121_Zscore <- LSCTE121_Zscore_calculation(opname="AMLbulk", path_counts_filtered="/results/AMLbulk/chromvar/", repeat_consensus_bed_path="/data/metadata/", out_dir=out_path)

# make heatmap annotation
annocol <- merge(metadata[colnames(zscore_AMLbulk),c("LSC17_score", "LSC_high_low","cyto")], as.data.frame(LSCTE121_Zscore[,"LSC_common_zscore"]), by=0) %>%
            column_to_rownames("Row.names")

colnames(annocol)[4] <- "LSCTE121_Z_score"
annocol$cyto <- gsub(" ", "", tolower(annocol$cyto))

# plot heatmap for LSC repeats
my.breaks <- c(seq(-10, 10, length.out=10))

cols <- colorRampPalette(brewer.pal(12, "Paired"))

mycolors <- cols(length(unique(annocol$cyto)))
names(mycolors) <- unique(annocol$cyto)
mycolors1 <- cols(length(unique(annocol$LSC_high_low)))
names(mycolors1) <- unique(annocol$LSC_high_low)

mycolors <- list(cyto=mycolors,LSC_high_low=mycolors1)

##----------------------
# Figure 4b
##----------------------
# plot for LSC repeats 
pdf(paste0(out_path, "Fig4b_Heatmap.AMLbulk.LSCrepeats_q0.01.pdf" ), onefile=F)
pheatmap(zscore_AMLbulk[LSC,],
breaks=my.breaks,
clustering_method="ward.D2",
clustering_distance_rows="euclidean",
clustering_distance_cols="euclidean",
col = brewer_yes,
annotation_col=annocol,
annotation_colors=mycolors,
cluster_rows=T, cluster_cols=T,
fontsize_col=5,fontsize_row=3,
fontsize =5,
show_rownames=F,
show_colnames=F,
scale="none")
dev.off()

# save metadata with LSCTE121 Z-score
meta <- merge(metadata, LSCTE121_Zscore, by=0)
write.csv(meta, paste0(out_path, "metadata_bulkAML_wLSCTE121.csv"), row.names=F)

##----------------------
# Suppl Figure 6a
##----------------------
# calculate LSCTE121 in LSC fractions
LSCTE121_Zscore <- LSCTE121_Zscore_calculation(opname="lsc", path_counts_filtered="/results/lsc/chromvar/", repeat_consensus_bed_path="/data/metadata/", out_dir=out_path)

data <- as.data.frame(merge(metadata, LSCTE121_Zscore, by=0))
data$eng <- ifelse(is.na(data$eng), "NE", data$eng)
data$LSC_Class_NE <- factor(data$eng, levels=c("high", "med", "low", "NE") )

ggplot(data, aes(x= LSC_Class_NE, y=LSC_common_zscore))+
        geom_boxplot() + theme_classic() + ylab("LSCTE121 Z-score") + xlab("LSC engraftment class") +
        geom_hline(yintercept = 0, color="grey", linetype="dotted")
ggsave("SupplFig6a_LSCTE121_Zscore_engraftment.pdf", path=out_path, device="pdf", height=5, width=7)

##----------------------
# Suppl Figure 6b
##----------------------
meta <- as.data.frame(meta)
ggscatter(meta, x = "LSC_common_zscore", y = "Blast_count_PB",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray"))+
        stat_cor(method = "pearson", label.x = -20, label.y = 250)  + xlab("LSCTE121 Z-score") + ylab("Blast count PB (x10e9/L) ") 
ggsave("SupplFig6b_LSCTE121_Zscore_Blast_count_correlation.pdf", path=out_path, device="pdf", height=5, width=7)

##----------------------
# Suppl Figure 6d
##----------------------

ggscatter(meta, x ="LSC17_score",  y = "LSC_common_zscore",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray"))+
        stat_cor(method = "pearson", label.y = -20, label.x = 1)  + ylab("LSCTE121 Z-score") + xlab("LSC17 score") 
ggsave("SupplFig6d_LSCTE121_Zscore_LSC17_score_correlation.pdf", path=out_path, device="pdf", height=5, width=7)


# ##----------------------
# # Suppl Figure 6f
# ##----------------------
# meta$cyto <- gsub(" ","", tolower(meta$cyto))
# ggplot(meta, aes(x=cyto, y=LSC_common_zscore, fill=cyto))+
#         geom_boxplot() + theme_classic() + xlab("Cytogenetics risk") + ylab("LSCTE121 Z-score") +
#         geom_hline(yintercept = 0, color="grey", linetype="dotted") +
#         scale_fill_manual(values=c("#F1C12B", "#C05E57", "#593C54"), na.value= "#6D6E71")
# ggsave("SupplFig6f_LSCTE121_Zscore_by_cytogentic_risk_groups.pdf", path=out_path, device="pdf", height=5, width=7)

