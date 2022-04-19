#----------------------------------------------------------------
## LSC analysis and Hemat comparisons
#----------------------------------------------------------------
## load dependencies
#----------------------------------------------------------------
rm(list = ls())
suppressMessages(library(chromVAR))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(dplyr))
suppressMessages(library(UpSetR))
suppressMessages(library(tibble))


source("/code/helper_scripts/chromvar.helperfunctions.R")
source("/code/helper_scripts/chromvar_analyse_functions.R")

out_path <- "/results/manuscript_figures/LSC/"
if(!dir.exists(out_path)){dir.create(out_path, recursive = TRUE)}

brewer_yes = c('#053061', '#2971B1', '#6AACD0', '#C1DDEB', '#F7F7F7', '#FACDB5', '#E58267', '#BB2933', '#67001F')

'%!in%' <- function(a,b){ !(a %in% b) }
#----------------------------------------------------------------
## Read in metadata
#----------------------------------------------------------------
metadata <- read.csv("/data/metadata/curated_metadata_lsc_hemat_new.csv", row.names=1)
metadata <- metadata[,-2]
metadata$CD34CD38 <- ifelse(metadata$CD34==1 & metadata$CD38==0, "CD34pos.CD38neg",
                        ifelse(metadata$CD34==1 & metadata$CD38==1, "CD34pos.CD38pos",
                        ifelse(metadata$CD34==0 & metadata$CD38==1, "CD34neg.CD38pos",
                        ifelse(metadata$CD34==0 & metadata$CD38==0, "CD34neg.CD38neg",NA))))


#------------------------------------------------------------------------------------
## LSC fractions and hematopoietic system -- LSC pos, neg, bulk and hematopoietic
#------------------------------------------------------------------------------------

##----------------------------
## Fig 2a
##----------------------------
zscore_lsc_hemat <- read.table("/results/lsc_hemat/chromvar/lsc_hemat.Zscore.txt",header=T, sep="\t", stringsAsFactors=F)
colnames(zscore_lsc_hemat) <- gsub("X", "", colnames(zscore_lsc_hemat))

annocol_lsc_hemat_15m <- data.frame(stringsAsFactors=F, group=metadata[colnames(zscore_lsc_hemat),c("V2", "V3")])
annocol_lsc_hemat_15m <- as.data.frame(apply(annocol_lsc_hemat_15m[1:2], 2, function(e) as.factor(e)))
annocol_lsc_hemat_15m$group.V2 <-factor(annocol_lsc_hemat_15m$group.V2, levels=c("LT-HSC","ST-HSC", "MLP", "CMP", "GMP", "MEP", "T Cell", "B Cell","Erythroid Precustor", "Monocyte", "Granulocyte", "Dendritic Cell", "NK Cell", "LSC.positive", "LSC.negative"))

my.breaks <- c(seq(-10, 10, length.out=10))
annocol1 <- annocol_lsc_hemat_15m

cols <- colorRampPalette(brewer.pal(12, "Paired"))
mycolors <- c(cols(13), "#D92727", "#0A3542")
names(mycolors) <- levels(annocol_lsc_hemat_15m$group.V2)
mycolors2 <- c("#F36F60", "#F79C76", "#A7D1D5","#D92727", "#0A3542")
names(mycolors2) <- unique(annocol_lsc_hemat_15m$group.V3)


mycolors <- list(group.V2=mycolors,group.V3= mycolors2)

pdf(paste0(out_path, "Fig2a_Heatmap.zscore_lsc_hemat_15m_TEonly.correlationclustering.pdf" ), onefile=F)
pheatmap(na.omit(zscore_lsc_hemat),
clustering_method="ward.D2",
breaks=my.breaks,
clustering_distance_rows="correlation",
clustering_distance_cols="correlation",
col = brewer_yes,
annotation_col=annocol1,
annotation_colors=mycolors,
cluster_rows=T, cluster_cols=T,
fontsize_col=5,fontsize_row=3,
fontsize =5,
show_rownames=F,
show_colnames=F,
scale="none")
dev.off()



#----------------------------------------------------------------
## lsc
#----------------------------------------------------------------
##----------------------------
## Fig 2b
##----------------------------
zscore_lsc <- read.table("/results/lsc/chromvar/lsc.Zscore.txt",header=T, sep="\t", stringsAsFactors=F)
colnames(zscore_lsc) <- gsub("X", "", colnames(zscore_lsc))


annocol_lsc <- data.frame(stringsAsFactors=F, group=metadata[colnames(zscore_lsc),])
annocol_lsc <- as.data.frame(apply(annocol_lsc[,c("group.V2","group.CD34CD38")], 2, function(e) as.factor(e)))

analysefunction(zscoresdf=zscore_lsc,
                dir="lsc",qval=0.01,
                annocol1=annocol_lsc,
                grouping=as.factor(annocol_lsc$group.V2),
                opname="lscpos_lscneg",
                out_path=out_path,
                fig="2b")

analysechromvar_medians(zscoresdf=zscore_lsc,
                dir="lsc",
                grouping=as.factor(annocol_lsc$group.V2),
                opname="lscpos_lscneg")




##----------------------------
## Suppl Fig 3a
##----------------------------
# plot Volcano plot for comparison 

lsc.diffdeviations <- read.table(paste0(out_path, "lsc.lscpos_lscneg.DiffDeviations.with.medians.txt"), header=T, sep="\t", stringsAsFactors=F)
lsc.diffdeviations$mediandiff <- lsc.diffdeviations$LSC.positive - lsc.diffdeviations$LSC.negative
lsc.diffdeviations$neglog10qval <- -log10(lsc.diffdeviations$p_value_adjusted)
lsc.diffdeviations$type <- ifelse(lsc.diffdeviations$p_value_adjusted < 0.01, "sig", "nonsig")
lsc.diffdeviations$type2 <- ifelse(lsc.diffdeviations$type=="sig" & lsc.diffdeviations$mediandiff>0, "LSC+",
                                                 ifelse(lsc.diffdeviations$type=="sig" & lsc.diffdeviations$mediandiff<0,"LSC-", "nonsig"))

pdf(paste0(out_path,"SupplFig3a_lsc.diffdeviations.volcanoplot.q0.01.pdf"), useDingbats=F)
ggplot(data=lsc.diffdeviations, aes(x=mediandiff, y=neglog10qval)) +
  geom_point(aes(color=type2)) +  
  #geom_text_repel(aes(label = epdf$Freqlabel)) +
  theme(text = element_text(size=12)) + 
  labs(x = "Median difference (LSC+ - LSC-)" , y="-log10(q-value)") + 
  scale_colour_manual(values=c(`LSC+`="#CC0033",nonsig="#000000", `LSC-`="#1e71aa")) +
  geom_hline(yintercept=2,linetype="dashed", size=0.5, col="grey") +
  xlim(-11, 11) + theme_classic()
dev.off()





#---------------------------
# Fig 2c
# family overview
#---------------------------
# load TE family annotation
exclude <- read.csv("/data/metadata/Repeats_nonTE.csv", stringsAsFactors=F)
repeat_metadata <- read.table("/data/metadata/repeat_metadata.txt", header=T, sep="\t", stringsAsFactors=F)
repfam <- read.table("/data/metadata/repeat_fam_mapping", header=F , sep="\t", stringsAsFactors=F)
repeat_metadata <- merge(repeat_metadata, repfam[,1:2], by.x="subfam",by.y="V2", all.x=T)
repeat_metadata[is.na(repeat_metadata)] <- "UNKNOWN"
repeat_metadata1 <- subset(repeat_metadata, repeat_metadata$repname %!in% exclude$repname)

# merge annotation with TE family annotation
lsc.diffdeviations$lsc_sig <- ifelse(lsc.diffdeviations$type == "sig", 1,0)
lsc.diffdeviations$lsc_sigstem <- ifelse(lsc.diffdeviations$type2 == "LSC+", 1,0)
lsc.diffdeviations$lsc_signonstem <- ifelse(lsc.diffdeviations$type2 == "LSC-", 1,0)

repfam_sig_lsc <- merge(lsc.diffdeviations, repeat_metadata, by="repname", all.x=TRUE, all.y=FALSE)



# make summary table to determine which families do not show significant subfamilies
summary_families_lsc <- plyr::ddply(repfam_sig_lsc, "V1", summarise, count=length(type), sig=sum(lsc_sig), lsc.pos =sum(lsc_sigstem), lsc.neg=sum(lsc_signonstem), .drop=FALSE)
write.csv(summary_families_lsc, paste0(out_path,"/LSCposLSCneg_repfam_overview.csv"))

# don't plot families without significant subfamilies
not <- summary_families_lsc[summary_families_lsc$sig == 0,][,1]

repfam_sig_lsc <- repfam_sig_lsc[!(repfam_sig_lsc$V1 %in% not),]
repfam_sig_lsc$famSort <- factor(repfam_sig_lsc$V1, )


pdf(paste0(out_path,"/Fig2c_lsc_pos_neg_repfam_overview_jitter_type2_onlysig_q0.01_TEonly.pdf" ), onefile=F, useDingbats=FALSE)
ggplot(repfam_sig_lsc,aes(x=famSort,y=mediandiff))+
    geom_jitter(aes(colour=type2), width=0.1)+
    coord_flip()+
    theme_bw() + 
    scale_x_discrete(limits = rev(levels(repfam_sig_lsc$famSort))) + 
    scale_colour_manual(values=c(`LSC+`="#CC0033",nonsig="#000000", `LSC-`="#1e71aa"))+ 
    labs(x = "TE family" , y="Median difference (Primitive - Mature)") 
dev.off()




##-------------------
# Fig 3d
# combined stacked bar plots
##-------------------
# load Hemat data
Hemat_stemprog_vs_diff.diffdeviations <- read.table(paste0("/results/manuscript_figures/Hemat/Hemat_TEonly.Hemat_stemprog_vs_diff.DiffDeviations.with.medians.txt"), header=T, sep="\t", stringsAsFactors=F)

Hemat_stemprog_vs_diff.diffdeviations$hemat_mediandiff <- Hemat_stemprog_vs_diff.diffdeviations$Primitive - Hemat_stemprog_vs_diff.diffdeviations$Mature
Hemat_stemprog_vs_diff.diffdeviations$hemat_neglog10qval <- -log10(Hemat_stemprog_vs_diff.diffdeviations$p_value_adjusted)
Hemat_stemprog_vs_diff.diffdeviations$hemat_type <- ifelse(Hemat_stemprog_vs_diff.diffdeviations$p_value_adjusted < 0.01, "sig", "nonsig")
Hemat_stemprog_vs_diff.diffdeviations$hemat_type2 <- ifelse(Hemat_stemprog_vs_diff.diffdeviations$hemat_type=="sig" & Hemat_stemprog_vs_diff.diffdeviations$hemat_mediandiff>0, "Primitive",
                                                 ifelse(Hemat_stemprog_vs_diff.diffdeviations$hemat_type=="sig" & Hemat_stemprog_vs_diff.diffdeviations$hemat_mediandiff<0,"Mature", "nonsig"))
Hemat_stemprog_vs_diff.diffdeviations.rep <- subset(Hemat_stemprog_vs_diff.diffdeviations, Hemat_stemprog_vs_diff.diffdeviations$p_value_adjusted < 0.01)[,5]

# make another df for the venn diagram/ percentage family plots comparison:
Hemat_stemprog_vs_diff.diffdeviations$hemat_sig <- ifelse(Hemat_stemprog_vs_diff.diffdeviations$hemat_type == "sig", 1,0)
Hemat_stemprog_vs_diff.diffdeviations$hemat_sigstem <- ifelse(Hemat_stemprog_vs_diff.diffdeviations$hemat_type2 == "Primitive", 1,0)
Hemat_stemprog_vs_diff.diffdeviations$hemat_signonstem <- ifelse(Hemat_stemprog_vs_diff.diffdeviations$hemat_type2 == "Mature", 1,0)


repfam_sig <- merge(Hemat_stemprog_vs_diff.diffdeviations[,c(5,8:12)], lsc.diffdeviations[,c(5,8:12)], by="repname", all=TRUE)

repfam_sig <- merge(repfam_sig, repeat_metadata1, by="repname", all=TRUE)
repfam_sig <- tidyr::replace_na(repfam_sig, list(hemat_type="nonsig", hemat_type2="nonsig", hemat_sig=0, hemat_sigstem=0,hemat_signonstem=0,type="nonsig",type2="nonsig", lsc_sig=0, lsc_sigstem=0, lsc_signonstem=0))

summary_families <- plyr::ddply(repfam_sig, "V1", summarise, Genomic=length(hemat_type), Primitive =sum(hemat_sigstem), Mature=sum(hemat_signonstem), `LSC+` =sum(lsc_sigstem), `LSC-`=sum(lsc_signonstem), .drop=FALSE)
fam_summary <- reshape::melt(summary_families, id.vars="V1")

#Family percentages stacked barplot
colourCount = length(unique(fam_summary$V1))
getPalette = colorRampPalette(brewer.pal(10, "RdYlBu"))
cols <- getPalette(colourCount)
names(cols) <- levels(factor(fam_summary$V1))

pdf(paste0(out_path,"Fig3d_family_count_plot_comparison.pdf"), onefile=F, useDingbats=FALSE)
ggplot(fam_summary, aes(fill=V1, y=value, x=variable)) + 
    geom_bar( stat="identity", position="fill") +coord_flip() +
    scale_fill_manual(values = getPalette(colourCount)) + theme_classic() + theme(legend.position="bottom", legend.box = "horizontal") +
    scale_x_discrete(limits = rev(levels(fam_summary$variable)))
dev.off() 





##-------------------
## Figure 2e
##-------------------
# Save significantly enriched TEs for Primitive, Mature, LSC+, LSC-
Primitive <-data.frame(repname = subset(repfam_sig, repfam_sig$hemat_sigstem ==1)[,1], group=rep("Primitive"))
Mature <- data.frame(repname = subset(repfam_sig, repfam_sig$hemat_signonstem ==1)[,1], group=rep("Mature"))
LSC_pos <-data.frame(repname = subset(repfam_sig, repfam_sig$lsc_sigstem ==1)[,1], group=rep("LSC+"))
LSC_neg <- data.frame(repname = subset(repfam_sig, repfam_sig$lsc_signonstem ==1)[,1], group=rep("LSC-"))
VennTable <- rbind(Primitive, Mature, LSC_pos, LSC_neg)
write.table(VennTable, paste0(out_path,"VennTable_TEonly.txt"))

pdf(paste0(out_path, "Fig2e_part1_UpsetR_plot.pdf"), onefile=F)
upset(fromList(list(Primitive=Primitive$repname, Mature=Mature$repname, `LSC+`=LSC_pos$repname, `LSC-`=LSC_neg$repname))) 
dev.off()

# part two
repfam_sig$comparison <- ifelse(repfam_sig$hemat_type2=="Mature" & repfam_sig$type2=="LSC-", "Mature_LSC-",
                                ifelse(repfam_sig$hemat_type2=="Primitive" & repfam_sig$type2=="LSC+", "Primitive_LSC+",
                                ifelse(repfam_sig$hemat_type2=="Primitive", "Primitive",
                                ifelse(repfam_sig$hemat_type2=="Mature", "Mature",
                                ifelse(repfam_sig$type2=="LSC+", "LSC+",
                                ifelse(repfam_sig$type2=="LSC-", "LSC-","nonsig"))))))

plot_comparison <- subset(repfam_sig, repfam_sig$comparison != "nonsig")
plot_comparison$comparison <- factor(plot_comparison$comparison, levels=c("Mature", "Primitive", "Mature_LSC-", "Primitive_LSC+","LSC-", "LSC+"))

pdf(paste0(out_path, "Fig2e_part2_barplot_family.pdf"), onefile=F, width=5, height=7)
ggplot(plot_comparison, aes(x=comparison, fill=V1)) + scale_fill_manual(values = cols)+
  stat_count() + theme_classic() + theme(legend.position="bottom", axis.text.x = element_text(angle = 90))
dev.off()



#---------------------------------
# Supplementary Fig 3b, c
# Cleaveland dotplots
#---------------------------------
# calculate cluster medians
zscore <- t(zscore_lsc_hemat) %>% data.frame()
colnames(zscore) <- gsub("^X", "", colnames(zscore))

metadata$anno <- ifelse(metadata$V3 %in% c("Hemat_stem", "Hemat_prog"), "Primitive", 
                                ifelse(metadata$V3 =="Hemat_diff", "Mature", 
                                ifelse(metadata$V3 =="CSC", "LSC+", "LSC-")))

zscore <- merge(zscore, metadata[,c("names","anno")],by.x=0, by.y="names", all.x = TRUE, all.y = FALSE)

zscore_grouped <- zscore[,2:973] %>% 
  group_by(anno) %>%
  summarise_all("median") %>%
  remove_rownames() %>%
  column_to_rownames("anno") %>%
  t() %>%
  as.data.frame()

# merge with repeat comparison annotation
zscore_plotting <- merge(zscore_grouped, plot_comparison[,c("repname","subfam","V1", "comparison")], by.x=0, by.y="repname", all.x = F) 

#PRIMITIVE - LSC+
# subset for Primitve - LSC+ plotting
primitive_LSCpos_zscore <- zscore_plotting %>% 
  select(c("Row.names","LSC+", "Primitive", "subfam", "V1", "comparison")) %>% 
  filter(comparison %in% c("Primitive", "Primitive_LSC+", "LSC+"))

# melt
primitive_LSCpos_zscore_plot <- reshape::melt(primitive_LSCpos_zscore, id.vars=c("Row.names","subfam","V1","comparison"))
colnames(primitive_LSCpos_zscore_plot)[1] <- "repname"


# plot the groups separately
primitive <- primitive_LSCpos_zscore_plot %>% filter(comparison %in% "Primitive")
LSCpos <- primitive_LSCpos_zscore_plot %>% filter(comparison %in% "LSC+")
LSCpos_primitive <- primitive_LSCpos_zscore_plot %>% filter(comparison %in% "Primitive_LSC+")

p1 <- ggplot(primitive, aes(value, subfam)) +
  geom_line(aes(group = subfam)) +
  geom_point(aes(color = variable)) +
  facet_grid(V1~., scales = "free_y", space= "free_y")  +
  theme_minimal() + scale_color_manual(values=c("#fcba03", "#fc5a03")) +
  theme(axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.text.y = element_text(angle = 0),
        legend.position="bottom") +
  ggtitle("Primitive repeats")

p2 <- ggplot(LSCpos, aes(value, subfam)) +
  geom_line(aes(group = subfam)) +
  geom_point(aes(color = variable)) +
  facet_grid(V1~., scales = "free_y", space= "free_y")  +
  theme_minimal() + scale_color_manual(values=c("#fcba03", "#fc5a03")) +
  theme(axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.text.y = element_text(angle = 0),
        legend.position="bottom") +
  ggtitle("LSC+ repeats")

p3 <- ggplot(LSCpos_primitive, aes(value, subfam)) +
  geom_line(aes(group = subfam)) +
  geom_point(aes(color = variable)) +
  facet_grid(V1~., scales = "free_y", space= "free_y")  +
  theme_minimal() + scale_color_manual(values=c("#fcba03", "#fc5a03")) +
  theme(axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.text.y = element_text(angle = 0),
        legend.position="bottom") +
  ggtitle("Shared LSC+/primitive repeats")

# combine plots
lay <- rbind(c(1,2),
             c(1,3))
ggsave("Suppl3bc_Cleaveland_LSCpos_primitve_repeats.pdf", plot=grid.arrange(grobs=list(p1,p2,p3), layout_matrix=lay), path=out_path, width=8,height=10, device="pdf", useDingbats=FALSE)



## MATURE - LSC-
# subset for Primitve - LSC+ plotting
mature_LSCneg_zscore <- zscore_plotting %>% 
  select(c("Row.names","LSC-", "Mature", "subfam", "V1", "comparison")) %>% 
  filter(comparison %in% c("Mature", "Mature_LSC-", "LSC-"))

# melt
mature_LSCneg_zscore_plot <- reshape::melt(mature_LSCneg_zscore, id.vars=c("Row.names","subfam","V1","comparison"))
colnames(mature_LSCneg_zscore_plot)[1] <- "repname"


mature <- mature_LSCneg_zscore_plot %>% filter(comparison %in% "Mature")
LSCneg <- mature_LSCneg_zscore_plot %>% filter(comparison %in% "LSC-")
LSCneg_MATURE <- mature_LSCneg_zscore_plot %>% filter(comparison %in% "Mature_LSC-")


p4 <- ggplot(mature, aes(value, subfam)) +
  geom_line(aes(group = subfam)) +
  geom_point(aes(color = variable)) +
  facet_grid(V1~., scales = "free_y", space= "free_y")  +
  theme_minimal() + scale_color_manual(values=c("#03bafc", "#0303fc")) +
  theme(axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.text.y = element_text(angle = 0),
        legend.position="bottom") +
  ggtitle("Mature repeats")

p5 <- ggplot(LSCneg, aes(value, subfam)) +
  geom_line(aes(group = subfam)) +
  geom_point(aes(color = variable)) +
  facet_grid(V1~., scales = "free_y", space= "free_y")  +
  theme_minimal() + scale_color_manual(values=c("#03bafc", "#0303fc")) +
  theme(axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.text.y = element_text(angle = 0),
        legend.position="bottom") +
  ggtitle("LSC- repeats")

p6 <- ggplot(LSCneg_MATURE, aes(value, subfam)) +
  geom_line(aes(group = subfam)) +
  geom_point(aes(color = variable)) +
  facet_grid(V1~., scales = "free_y", space= "free_y")  +
  theme_minimal() + scale_color_manual(values=c("#03bafc", "#0303fc")) +
  theme(axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.text.y = element_text(angle = 0),
        legend.position="bottom") +
  ggtitle("Shared LSC-/mature repeats")


# combine plots
lay <- rbind(c(1,2),
             c(1,3),
             c(1,3))
ggsave("Suppl3bc_Cleaveland_LSCneg_mature_repeats.pdf", plot=grid.arrange(grobs=list(p4,p5,p6), layout_matrix=lay), path=out_path, width=8,height=12, device="pdf", useDingbats=FALSE)
