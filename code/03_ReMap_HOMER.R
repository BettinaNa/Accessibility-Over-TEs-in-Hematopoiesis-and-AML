## LOLA and HOMER analysis 

rm(list = ls())

#library(Matrix)
suppressMessages(library(ggrepel))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(UpSetR))
suppressMessages(library(magick))
suppressMessages(library(ggpubr))
suppressMessages(library(grid))


out_path <- "/results/manuscript_figures/LOLA_HOMER/"
dir.create(out_path)

## TOP 2.5% TFs

comparisons<- c("Primitive.only","LSCpos.only","Primitive.LSCpos")
figs <- c("SupplFig4a","Fig3a","SupplFig4b")
top <- list()
RUNX1 <- list()
ERG <- list() 


for (i in seq_along(comparisons)){
  significant_TFs <- read.table(paste0("/results/LOLA/",comparisons[i],"/combined_lola.txt"), header=T, sep="\t", stringsAsFactors=F)
  # calculate number of times each TF was enriched
  top[[eval(comparisons[i])]] <- significant_TFs %>% 
            count(antibody) %>% 
            arrange(-n) %>% 
            mutate(id = row_number(), top0.025 = ifelse(max(id)*0.025>=id, "top2.5percent", "nr"), min_count=min(n[which(top0.025=="top2.5percent")])) %>% 
            filter(n >= min_count) %>%
            select(antibody, n) %>%
            rename(freq=n)

  # plot 
  pdf(paste0(out_path, figs[i], "_", comparisons[i], "_freq.top2.5pc.graph.pdf"), onefile=F, height=5, width=4)
    print(ggplot(top[[eval(comparisons[i])]],aes(x=reorder(antibody, freq),y=freq))+
    geom_point(stat="identity", size = 5, stroke = 0, shape=16)+
    theme_bw()+
    coord_flip() +
    labs(x="Transcription Factors", y="Number of TE subfamilies") +
    ggtitle(paste0("Top 2.5% Transcription Factors at ",comparisons[i]," TE subfamilies")) +
    theme(axis.text = element_text(size = 12)))
  dev.off()

  
  RUNX1[[eval(comparisons[i])]] <- significant_TFs %>% 
                                    filter(antibody == "RUNX1") %>%
                                    select(antibody, TE_subfamily)
  
  ERG[[eval(comparisons[i])]] <- significant_TFs %>% 
                                    filter(antibody == "ERG") %>%
                                    select(antibody, TE_subfamily)

}


# TF comparison Primitive, LSCpos, both
#-------------------
# Figure 3b
#-------------------
pdf(paste0(out_path, "Fig3b_UpsetR.Comparing_Primitive.LSCp.both.TFs.top2.5pc.graph.pdf"), useDingbats=F,onefile=F)
upset(fromList(list(Primitive=top[["Primitive.only"]]$antibody, `LSC positive`=top[["LSCpos.only"]]$antibody, `Primitive & LSC positive`=top[["Primitive.LSCpos"]]$antibody)),order.by = "freq")
dev.off()


#-------------------
# Suppl Figure 4c
#-------------------
pdf(paste0(out_path, "SupplFig4c_part1_UpsetR.Primitive_ERG_RUNX1_subfamilies.pdf"), useDingbats=F,onefile=F)
upset(fromList(list(`ERG enriched`=ERG[["Primitive.only"]]$TE_subfamily, `RUNX1 enriched`=RUNX1[["Primitive.only"]]$TE_subfamily)),order.by = "freq")
dev.off()

pdf(paste0(out_path, "SupplFig4c_part2_UpsetR.LSCpos.only_ERG_RUNX1_subfamilies.pdf"), useDingbats=F,onefile=F)
upset(fromList(list(`ERG enriched`=ERG[["LSCpos.only"]]$TE_subfamily, `RUNX1 enriched`=RUNX1[["LSCpos.only"]]$TE_subfamily)),order.by = "freq")
dev.off()

pdf(paste0(out_path, "SupplFig4c_part3_UpsetR.Primitive.LSCpos_ERG_RUNX1_subfamilies.pdf"), useDingbats=F,onefile=F)
upset(fromList(list(`ERG enriched`=ERG[["Primitive.LSCpos"]]$TE_subfamily, `RUNX1 enriched`=RUNX1[["Primitive.LSCpos"]]$TE_subfamily)),order.by = "freq")
dev.off()


#------------------
# Figure 1e
#------------------

My_data <- read.table("/results/LOLA/LT.HSC/1032_MER57B1_remap_2018/allEnrichments.tsv", header=T, check.names=F,sep="\t", stringsAsFactors=F)
My_data$neglog10qval <- -log10(My_data$qValue)
pdf(paste0(out_path,"Fig1e_part1_ALL_TFs.Remap_hg38_1032_MER57B1_remap_2018.pdf" ), onefile=F)
ggplot(My_data, aes(x = oddsRatio, y = neglog10qval, label=description)) +
geom_point(size=3,stroke = 0, shape=16) +
theme_bw() +
labs(title = "TFs on HERV1_MER57B1 (LT-HSCs)",x = "OddsRatio", y = "-log10(q-value)") +
geom_vline(xintercept = 1, linetype = "dashed", colour = "firebrick2") + geom_hline(yintercept = 1.3, linetype = "dashed", colour = "firebrick2") +
geom_text(aes(label=ifelse(neglog10qval>1.3,as.character(description),'')),hjust=0.5,vjust=-0.5)
dev.off()

My_data <- read.table("/results/LOLA/LT.HSC/1033_MER57B2_remap_2018/allEnrichments.tsv", header=T, check.names=F,sep="\t", stringsAsFactors=F)
My_data$neglog10qval <- -log10(My_data$qValue)
pdf(paste0(out_path,"Fig1e_part2_ALL_TFs.Remap_hg38_1033_MER57B2_remap_2018.pdf" ), onefile=F)
ggplot(My_data, aes(x = oddsRatio, y = neglog10qval, label=description)) +
geom_point(size=3,stroke = 0, shape=16) +
theme_bw() +
labs(title = "TFs on HERV1_MER57B2 (LT-HSCs)",x = "OddsRatio", y = "-log10(q-value)") +
geom_vline(xintercept = 1, linetype = "dashed", colour = "firebrick2") + geom_hline(yintercept = 1.3, linetype = "dashed", colour = "firebrick2") +
geom_text(aes(label=ifelse(neglog10qval>1.3,as.character(description),'')),hjust=0.5,vjust=-0.5)
dev.off()

My_data <- read.table("/results/LOLA/LT.HSC/826_LTR67B_remap_2018/allEnrichments.tsv", header=T, check.names=F,sep="\t", stringsAsFactors=F)
My_data$neglog10qval <- -log10(My_data$qValue)
pdf(paste0(out_path,"Fig1e_part3_ALL_TFs.Remap_hg38_826_LTR67B_remap_2018.pdf" ), onefile=F)
ggplot(My_data, aes(x = oddsRatio, y = neglog10qval, label=description)) +
geom_point(size=3,stroke = 0, shape=16) +
theme_bw() +
labs(title = "TFs on HERV3_LTR67B (LT-HSCs)",x = "OddsRatio", y = "-log10(q-value)") +
geom_vline(xintercept = 1, linetype = "dashed", colour = "firebrick2") + geom_hline(yintercept = 1.3, linetype = "dashed", colour = "firebrick2") +
geom_text(aes(label=ifelse(neglog10qval>5,as.character(description),'')),hjust=0.5,vjust=-0.5)
dev.off()

My_data <- read.table("/results/LOLA/LT.HSC/897_MamGypsy2-LTR_remap_2018/allEnrichments.tsv", header=T, check.names=F,sep="\t", stringsAsFactors=F)
My_data$neglog10qval <- -log10(My_data$qValue)
pdf(paste0(out_path,"Fig1e_part4_ALL_TFs.Remap_hg38_897_MamGypsy2-LTR_remap_2018.pdf" ), onefile=F)
ggplot(My_data, aes(x = oddsRatio, y = neglog10qval, label=description)) +
geom_point(size=3,stroke = 0, shape=16) +
theme_bw() +
labs(title = "TFs on LTR(no HERVs)_MamGypsy2-LTR (LT-HSCs)",x = "OddsRatio", y = "-log10(q-value)") +
geom_vline(xintercept = 1, linetype = "dashed", colour = "firebrick2") + geom_hline(yintercept = 1.3, linetype = "dashed", colour = "firebrick2") +
geom_text(aes(label=ifelse(neglog10qval>1.3,as.character(description),'')),hjust=0.5,vjust=-0.5)
dev.off()



#-----------------
# Suppl Fig 4d
#-----------------

# Primitive only

My_data <- read.table("/results/LOLA/Primitive.only/843_LTR78_remap_2018/allEnrichments.tsv", header=T, check.names=F,sep="\t", stringsAsFactors=F)
My_data$neglog10qval <- -log10(My_data$qValue)
pdf(paste0(out_path,"SupplFig4d_part1_ALL_TFs.Remap_hg38_843_LTR78_remap_2018.pdf" ), onefile=F)
ggplot(My_data, aes(x = oddsRatio, y = neglog10qval, label=description)) +
geom_point(size=3,stroke = 0, shape=16) +
theme_bw() +
labs(title = "TFs on HERV1_LTR78 (Primitive)",x = "OddsRatio", y = "-log10(q-value)") +
geom_vline(xintercept = 1, linetype = "dashed", colour = "firebrick2") + geom_hline(yintercept = 1.3, linetype = "dashed", colour = "firebrick2") +
geom_text(aes(label=ifelse(neglog10qval>=1.3,as.character(description),'')),hjust=0.5,vjust=-0.5)
dev.off()

My_data <- read.table("/results/LOLA/Primitive.only/782_LTR41_remap_2018/allEnrichments.tsv", header=T, check.names=F,sep="\t", stringsAsFactors=F)
My_data$neglog10qval <- -log10(My_data$qValue)
pdf(paste0(out_path,"SupplFig4d_part4_ALL_TFs.Remap_hg38_4782_LTR41_remap_2018.pdf" ), onefile=F)
ggplot(My_data, aes(x = oddsRatio, y = neglog10qval, label=description)) +
geom_point(size=3,stroke = 0, shape=16) +
theme_bw() +
labs(title = "TFs on HERV3_LTR41 (Primitive)",x = "OddsRatio", y = "-log10(q-value)") +
geom_vline(xintercept = 1, linetype = "dashed", colour = "firebrick2") + geom_hline(yintercept = 1.3, linetype = "dashed", colour = "firebrick2") +
geom_text(aes(label=ifelse(neglog10qval>=40,as.character(description),'')),hjust=0.5,vjust=-0.5)
dev.off()

# Primitive and LSC+

My_data <- read.table("/results/LOLA/Primitive.LSCpos/1157_MLT1E2_remap_2018/allEnrichments.tsv", header=T,check.names=F,sep="\t", stringsAsFactors=F)
My_data$neglog10qval <- -log10(My_data$qValue)
pdf(paste0(out_path,"SupplFig4d_part3_ALL_TFs.Remap_hg38_1157_MLT1E2_remap_2018.pdf" ), onefile=F)
ggplot(My_data, aes(x = oddsRatio, y = neglog10qval, label=description)) +
geom_point(size=3,stroke = 0, shape=16) +
theme_bw() +
labs(title = "TFs on HERV3_MLT1E2 (Primitive - LSC+)",x = "OddsRatio", y = "-log10(q-value)") +
geom_vline(xintercept = 1, linetype = "dashed", colour = "firebrick2") + geom_hline(yintercept = 1.3, linetype = "dashed", colour = "firebrick2") +
geom_text(aes(label=ifelse(neglog10qval>1.3,as.character(description),'')),hjust=0.5,vjust=-0.5)
dev.off()

# LSC+

My_data <- read.table("/results/LOLA/LSCpos.only/826_LTR67B_remap_2018/allEnrichments.tsv", header=T, check.names=F,sep="\t", stringsAsFactors=F)
My_data$neglog10qval <- -log10(My_data$qValue)
pdf(paste0(out_path,"SupplFig4d_part4_ALL_TFs.Remap_hg38_826_LTR67B_remap_2018.pdf" ), onefile=F)
ggplot(My_data, aes(x = oddsRatio, y = neglog10qval, label=description)) +
geom_point(size=3,stroke = 0, shape=16) +
theme_bw() +
labs(title = "TFs on HERV3_LTR67B (LSC+)",x = "OddsRatio", y = "-log10(q-value)") +
geom_vline(xintercept = 1, linetype = "dashed", colour = "firebrick2") + geom_hline(yintercept = 1.3, linetype = "dashed", colour = "firebrick2") +
geom_text(aes(label=ifelse(neglog10qval>=5.35,as.character(description),'')),hjust=0.5,vjust=-0.5)
dev.off()


#-----------------
# Figure 1f
#-----------------
LT.HSCs.specific.TEs.ALL.motif.5sf <- read.delim("/data/homer_data/LT.HSCs.specific.TEs.ALL.motif.5sf.txt")
dt<-LT.HSCs.specific.TEs.ALL.motif.5sf


positions <- c("MER57B1", "MER57B2", "LTR67B", "MamGypsy2-LTR")

a<-ggplot(data=dt[1:4,],aes(x=TE_subfamily,-(-log10(q.value..Benjamini.)),y=-log10(q.value..Benjamini.)))+
  geom_col(width=0.5)+scale_y_continuous(limits = c(0, NA))+scale_x_discrete(limits = positions)+theme_minimal()+ 
  labs(title ="CTCF(Zf)/CD4+-CTCF",x="TE Subfamily",y="-log10(q-value)")+theme(axis.text.x=element_text(angle=30, hjust=1))+ 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")


b<-ggplot(data=dt[5:8,],aes(x=TE_subfamily,-(-log10(q.value..Benjamini.)),y=-log10(q.value..Benjamini.)))+
  geom_col(width=0.5)+scale_y_continuous(limits = c(0, NA))+scale_x_discrete(limits = positions)+theme_minimal()+
  labs(title ="ERG(ETS)/VCaP-ERG",x="TE Subfamily",y="-log10(q-value)")+theme(axis.text.x=element_text(angle=30, hjust=1))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")


c<-ggplot(data=dt[9:12,],aes(x=TE_subfamily,-(-log10(q.value..Benjamini.)),y=-log10(q.value..Benjamini.)))+
  geom_col(width=0.5)+scale_y_continuous(limits = c(0, NA))+scale_x_discrete(limits = positions)+theme_minimal()+
  labs(title ="HOXA9/HSC-Hoxa9",x="TE Subfamily",y="-log10(q-value)")+theme(axis.text.x=element_text(angle=30, hjust=1))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")


d<-ggplot(data=dt[13:16,],aes(x=TE_subfamily,-(-log10(q.value..Benjamini.)),y=-log10(q.value..Benjamini.)))+
  geom_col(width=0.5)+scale_y_continuous(limits = c(0, NA))+scale_x_discrete(limits = positions)+theme_minimal()+
  labs(title ="RUNX1(Runt)/Jurkat-RUNX1",x="TE Subfamily",y="-log10(q-value)")+theme(axis.text.x=element_text(angle=30, hjust=1))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")


figure <- ggarrange(a,b,c,d,ncol = 2, nrow = 2)

pdf(paste0(out_path,"Fig1f_HOXA9_ERG_RUNX1_CTCF_LT.HSC.pdf" ), onefile=F)
print(figure)
dev.off()


#-----------------
# Figure 3c
#-----------------
#NFY
motifs<- read.delim("/data/homer_data/EXCLUSIVE.LSCp.TEs.NFYA.enriched.NFYA.motifs.txt", stringsAsFactors=F)
dt<-motifs
uniq_motifs<-unique(dt$Motif.Name)

positions <- sort(dt$TE_subfamily)

a2 <- ggplot(data=dt,aes(x=TE_subfamily,-(-log10(q.value..Benjamini.)),y=-log10(q.value..Benjamini.)))+geom_col(width=0.3)+scale_y_continuous(limits = c(0, NA))+scale_x_discrete(limits = positions)+theme_minimal()+ labs(title=uniq_motifs,x="TE Subfamily",y="-log10(q-value)")+theme(text = element_text(size=5),axis.text.x=element_text(angle=70, hjust=1))+ geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

#FLI1
motifs<- read.delim("/data/homer_data/EXCLUSIVE.LSCp.TEs.FLI1.enriched.FLI1.motifs.txt", stringsAsFactors=F)
dt<- motifs %>% filter(Motif.Name == "EWS:FLI1-fusion(ETS)/SK_N_MC-EWS:FLI1")

uniq_motifs<-unique(dt$Motif.Name)
positions <- sort(dt$TE_subfamily)

a1 <- ggplot(data=dt,aes(x=TE_subfamily,-(-log10(q.value..Benjamini.)),y=-log10(q.value..Benjamini.)))+geom_col(width=0.3)+scale_y_continuous(limits = c(0, NA))+scale_x_discrete(limits = positions)+theme_minimal()+ labs(title=uniq_motifs,x="TE Subfamily",y="-log10(q-value)")+theme(text = element_text(size=5),axis.text.x=element_text(angle=70, hjust=1))+ geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

figure <- ggarrange(a1,a2,ncol = 1, nrow = 2)

pdf(paste0(out_path,"Fig3c_FLI1_fusion_NFY_LSCpos.pdf" ), onefile=F)
print(figure)
dev.off()

#-----------------
# Suppl Figure 4e - ERG
#-----------------
motifs<- read.delim("/data/homer_data/ALL.TEs.ERG.enriched.ERG.motifs.txt", stringsAsFactors=F)
dt<- motifs %>% filter(Motif.Name =="ERG(ETS)/VCaP-ERG")

uniq_motifs<-unique(dt$Motif.Name)
positions <- sort(dt$TE_subfamily)

a1 <- ggplot(data=dt,aes(x=TE_subfamily,-(-log10(q.value..Benjamini.)),y=-log10(q.value..Benjamini.)))+geom_col(width=0.3)+scale_y_continuous(limits = c(0, NA))+scale_x_discrete(limits = positions)+theme_minimal()+ labs(title=uniq_motifs,x="TE Subfamily",y="-log10(q-value)")+theme(text = element_text(size=5),axis.text.x=element_text(angle=70, hjust=1))+ geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

figure <- ggarrange(a1,ncol = 1, nrow = 1)

pdf(paste0(out_path,"SupplFig4e_ERG_all_Primitive_LSC+.pdf" ), onefile=F, height=4, width=12)
print(figure)
dev.off()



#-----------------
# Suppl Figure 4f - RUNX1
#-----------------

motifs<- read.delim("/data/homer_data/ALL.TEs.RUNX1.enriched.RUNX1.motifs.txt", stringsAsFactors=F)
dt<-motifs %>% filter(Motif.Name =="RUNX1(Runt)/Jurkat-RUNX1")

uniq_motifs<-unique(dt$Motif.Name)

positions <- sort(dt$TE_subfamily)

a1 <- ggplot(data=dt,aes(x=TE_subfamily,-(-log10(q.value..Benjamini.)),y=-log10(q.value..Benjamini.)))+geom_col(width=0.3)+scale_y_continuous(limits = c(0, NA))+scale_x_discrete(limits = positions)+theme_minimal()+ labs(title=uniq_motifs,x="TE Subfamily",y="-log10(q-value)")+theme(text = element_text(size=5),axis.text.x=element_text(angle=70, hjust=1))+ geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

figure <- ggarrange(a1,ncol = 1, nrow = 1)

pdf(paste0(out_path,"SupplFig4f_RUNX1_all_Primitive_LSC+.pdf" ), onefile=F, height=4, width=12)
print(figure)
dev.off()


#-----------------
# Suppl Figure 4g - CTCF
#-----------------
motifs<- read.delim("/data/homer_data/Primitive.TES.CTCF.enriched.CTCF.motifs.txt", stringsAsFactors=F)
dt<-motifs %>% filter(Motif.Name == "CTCF(Zf)/CD4+-CTCF")

uniq_motifs<-unique(dt$Motif.Name)

positions <- sort(dt$TE_subfamily)

a1 <- ggplot(data=dt,aes(x=TE_subfamily,-(-log10(q.value..Benjamini.)),y=-log10(q.value..Benjamini.)))+geom_col(width=0.3)+scale_y_continuous(limits = c(0, NA))+scale_x_discrete(limits = positions)+theme_minimal()+ labs(title=uniq_motifs,x="TE Subfamily",y="-log10(q-value)")+theme(text = element_text(size=5),axis.text.x=element_text(angle=70, hjust=1))+ geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

figure <- ggarrange(a1,ncol = 1, nrow = 1)

pdf(paste0(out_path,"SupplFig4g_CTCF_Primitive.pdf" ), onefile=F, height=4, width=12)
print(figure)
dev.off()








