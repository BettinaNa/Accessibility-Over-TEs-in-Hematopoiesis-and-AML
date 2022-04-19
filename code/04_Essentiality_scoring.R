# Essentiality scores of LSC+ specific transcription regulators - DepMap - CRISPR and RNAi data
# complete data sets available from DepMap downloads
# CRISPR: https://ndownloader.figshare.com/files/26261293
# RNAi: https://ndownloader.figshare.com/files/13515395

rm(list = ls())
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
options(scipen = 999)

#----------------------------------------------
# Part 1: CRISPR
#----------------------------------------------
# define out_path
out_path <- "/results/manuscript_figures/Essentiality_score/CRISPR/"
if(!dir.exists(out_path)){dir.create(out_path, recursive = TRUE)}


# make loop to go through each of the genes
TFs <- c("FLI1", "GABPA", "LYL1", "MEIS1", "NFYA", "RUNX1T1", "TCF3")
bestID <- "AML"
combined.ALL <- data.frame()

for (t in TFs){
  data <- read.table(paste0("/data/essentiality_data/CRISPR/",eval(t),"_CRISPR_all_essentiality_scores.txt"), header=T, sep="\t", stringsAsFactors=F)
  data_5_ST <- data[data$Lineage.subtype %in% names(which(table(data$Lineage.subtype) > 4)), ] %>% as.data.frame()

  medians_data <- data_5_ST %>% 
                    group_by(Lineage.subtype) %>% 
                    summarise(Median=median(!!sym(eval(t)))) %>% 
                    as.data.frame()
 
  # Generate file with BH corrected p-value
  # pair-wise t-test with BH correction - STATS WERE ADDED MANUALLY ON THE GRAPH AFTERWARDS
  stat <- pairwise.t.test(data_5_ST[,eval(t)], data_5_ST$Lineage.subtype,
                  p.adjust.method = "BH")
  #
  stat_AML.row <- stat$p.value[as.character(bestID),] %>% as.data.frame() %>% tibble::rownames_to_column("Lineage.subtype") %>% rename(.,stat_AML=.) %>% na.omit()
  stat_AML.col <- stat$p.value[,as.character(bestID)] %>% as.data.frame() %>% tibble::rownames_to_column("Lineage.subtype") %>% rename(.,stat_AML=.) %>% na.omit()
  #
  stat.AML_ALL <- rbind(stat_AML.row, stat_AML.col)
  #
  combined.data <- merge(medians_data, stat.AML_ALL, by="Lineage.subtype", all.x=TRUE)
  combined.data$Gene <- c(eval(t))

  combined.ALL <- rbind(combined.ALL, combined.data)

}


combined.ALL$Lineage.subtype <- chartr("_", " ", combined.ALL$Lineage.subtype)

combined.ALL$stat_AML <- ifelse(combined.ALL$stat_AML > 0.05, "ns",
                            ifelse(combined.ALL$stat_AML < 0.05 & combined.ALL$stat_AML > 0.01, "0.05",
                            ifelse(combined.ALL$stat_AML < 0.01 & combined.ALL$stat_AML > 0.001, "0.01", "0.001")))

ordered <- levels(factor(combined.ALL$Lineage.subtype))[c(3,2,4:36)]
combined.ALL$cancer_type <- factor(combined.ALL$Lineage.subtype, levels =ordered)

combined.ALL$regulator <- factor(combined.ALL$Gene, levels=c("TCF3", "RUNX1T1", "MEIS1", "NFYA", "GABPA", "LYL1", "FLI1"))

pdf(paste0(out_path, "Figure3d_Grid.essentiality.scores.stat.ALL.cell.lines.CRISPR.pdf" ), onefile=F)
ggplot(combined.ALL, aes(x=cancer_type, y=regulator)) +
  geom_tile(aes(fill = Median, color=stat_AML, width=0.7, height=0.5), size=0.75) +
  theme_bw() +
  theme(legend.position="bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Cancer types" , y="Transcription regulators", fill="Median essentiality score", color = "adjusted p-value") +
  scale_color_manual(values=c('ns'="gray50", '0.05'="lightpink1", '0.01'="tomato", '0.001'="red")) +
  scale_fill_gradient(low="black", high="white")

dev.off()


##------------------
## Suppl. Fig. 5a
##------------------
# We extracted the data of all AML cell lines from the DepMap data base ("ACH-000004", "ACH-000005","ACH-000045","ACH-000113","ACH-000146","ACH-000168","ACH-000198","ACH-000263","ACH-000294",
# "ACH-000336","ACH-000362","ACH-000387","ACH-000406","ACH-000487","ACH-000498","ACH-000557","ACH-000602","ACH-000770","ACH-001129","ACH-001574","ACH-001577","ACH-001618","ACH-001647") to compare # essentiality scores of each of LSC+ transcription regulators in AML cell lines and merge the data in one file: Combined.ALL.GENES.AML.lines.only.txt

ALL_essentiality_AML <- read.table("/data/essentiality_data/CRISPR/Combined.ALL.GENES.AML.lines.only.txt", header=T, sep="\t", stringsAsFactors=F)
plotlist <- list()

colors <- c("grey50", "deepskyblue1")

for (t in TFs){
  
  ALL_essentiality_AML1 <- ALL_essentiality_AML
  ALL_essentiality_AML1$Group <- ifelse(str_detect(ALL_essentiality_AML1$GeneID, paste0(eval(t), "\\.")),eval(t),"All_Genes")
  
  sub <- subset(ALL_essentiality_AML1, ALL_essentiality_AML1$Group == eval(t))
  names(colors) <- c("All_Genes", eval(t))
  plotlist[[eval(t)]] <- ggplot(ALL_essentiality_AML1, aes(x=Group, y=Essentiality_score, fill=Group)) +
                          geom_violin(draw_quantiles = c(0.5)) +
                          theme(legend.position="none") +
                          geom_jitter(data=sub, aes(x=Group, y=Essentiality_score), color="black", shape=16, size=2, position=position_jitter(0.2)) +
                          labs(title=paste0("Essentiality scores - ",eval(t)), x = "Genes" , y="Essentiality scores") +
                          scale_fill_manual(values=colors)

  # write file with stats
  NoNAs <- na.omit(ALL_essentiality_AML1)
  stats <- broom::tidy(t.test(NoNAs$Essentiality_score, sub$Essentiality_score, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F))
  write.table(stats, paste0(out_path, "T.test_", eval(t),".txt"))

}

pdf(paste0(out_path, "SupplFig5a_Violin.essentiality.scores.AML.cell.lines.all.genes.CRISPR.pdf" ), onefile=F, height=8, width=16)
plot_grid(plotlist=plotlist, ncol=4, nrow=2)
dev.off()




#----------------------------------------------
# Part 2: RNAi
#----------------------------------------------

# define out_path
out_path <- "/results/manuscript_figures/Essentiality_score/RNAi/"
if(!dir.exists(out_path)){dir.create(out_path, recursive = TRUE)}


# make loop to go through each of the TFs
TFs <- c("FLI1", "LYL1", "MEIS1", "NFYA", "RUNX1T1", "TCF3")
bestID <- "AML"
combined.ALL <- data.frame()

ALL_essentiality_all_lines <- read.table("/data/essentiality_data/RNAi/Six_genes_all_essentiality_scores.v4.txt", header=T, sep="\t", stringsAsFactors=F)
ALL_essentiality_all_lines_5 <- ALL_essentiality_all_lines[ALL_essentiality_all_lines$Group %in% names(which(table(ALL_essentiality_all_lines$Group) > 4)), ] %>% as.data.frame()


for (t in TFs){
   data <- ALL_essentiality_all_lines_5 %>% 
                    select(Cell.line, !!sym(eval(t)),Group) %>% 
                    as.data.frame()

   medians_data <- data %>%
                    group_by(Group) %>% 
                    summarise(Median=median(!!sym(eval(t)))) %>% 
                    as.data.frame()

  # Generate file with BH corrected p-value
  # pair-wise t-test with BH correction - STATS WERE ADDED MANUALLY ON THE GRAPH AFTERWARDS
 
  stat <- pairwise.t.test(data[,eval(t)], data$Group,
                  p.adjust.method = "BH")
  #
  stat_AML.col <- stat$p.value[,as.character(bestID)] %>% as.data.frame() %>% tibble::rownames_to_column("Group") %>% rename(.,stat_AML=.) %>% na.omit()
  #
  combined.data <- merge(medians_data, stat_AML.col, by="Group", all.x=TRUE)
  combined.data$Gene <- c(eval(t))

  combined.ALL <- rbind(combined.ALL, combined.data)

}


combined.ALL$stat_AML <- ifelse(combined.ALL$stat_AML > 0.05, "ns",
                            ifelse(combined.ALL$stat_AML < 0.05 & combined.ALL$stat_AML > 0.01, "0.05",
                            ifelse(combined.ALL$stat_AML < 0.01 & combined.ALL$stat_AML > 0.001, "0.01", "0.001")))

combined.ALL$regulator <- factor(combined.ALL$Gene, levels=c("TCF3", "RUNX1T1", "MEIS1", "NFYA", "LYL1", "FLI1"))

##--------------
# Figure 3e
##--------------

pdf(paste0(out_path,"Fig3e_Grid.essentiality.scores.stat.ALL.cell.lines.RNAi.pdf" ), onefile=F)
ggplot(combined.ALL, aes(x=Group, y=regulator)) +
  geom_tile(aes(fill = Median, color=stat_AML, width=0.7, height=0.5), size=0.75) +
  theme_bw() +
  theme(legend.position="bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Cancer types" , y="Transcription regulators", fill="Median essentiality score", color = "adjusted p-value") +
  scale_color_manual(values=c('ns'="gray50", '0.05'="lightpink1", '0.01'="tomato", '0.001'="red")) +
  scale_fill_gradient(low="black", high="white")
dev.off()


##--------------
# Suppl. Fig. 5b
##--------------

plotlist <- list()
colors <- c("grey50", "deepskyblue1")

for (t in TFs){
  ALL_essentiality_AML <- read.table(paste0("/data/essentiality_data/RNAi/ALL.essentiality.scores.AML.lines.",eval(t),".txt"), header=T, sep="\t", stringsAsFactors=F)
  sub <- subset(ALL_essentiality_AML, ALL_essentiality_AML$Group == eval(t))

  names(colors) <- c("ALL_GENES", eval(t))
  plotlist[[eval(t)]] <- ggplot(ALL_essentiality_AML, aes(x=Group, y=AML_cell_line, fill=Group)) +
                          geom_violin(draw_quantiles = c(0.5)) +
                          theme(legend.position="none") +
                          geom_jitter(data=sub, aes(x=Group, y=AML_cell_line), color="black", shape=16, size=2, position=position_jitter(0.2)) +
                          labs(title=paste0("Essentiality scores - ",eval(t)), x = "Genes" , y="Essentiality scores") +
                          scale_fill_manual(values=colors)
    
  # write file with stats
  NoNAs <- na.omit(ALL_essentiality_AML)
  stats <- broom::tidy(t.test(NoNAs$AML_cell_line, sub$AML_cell_line, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F))
  write.table(stats, paste0(out_path, "T.test_", eval(t),".txt"))

}

pdf(paste0(out_path, "SupplFig5b_Violin.essentiality.scores.AML.cell.lines.all.genes.RNAi.pdf" ), onefile=F, height=8, width=12)
plot_grid(plotlist=plotlist, ncol=3, nrow=2)
dev.off()


