##-------------------------------------------------------------
## AML survival analysis 
## using different clinical survival endpoints - overall survival (OS), Disease-Free_interval (DFI)
##-------------------------------------------------------------
rm(list = ls())
suppressMessages(library(survival))

source("/code/helper_scripts/KM.survival.R")
source("/code/helper_scripts/COX.regression_multi.R")
run_km <- function(DATASET,OPNAME,COLN) { KM.survival(x=1,DATA=DATASET,y=OPNAME,NAME=COLN) }

out_path <- "/results/manuscript_figures/AMLbulk/"
opname<-"zscore_AMLbulk"

metadata_groups <- read.csv(paste0(out_path,"metadata_bulkAML_wLSCTE121.csv"),row.names=1, stringsAsFactors=F)


if(!dir.exists(paste0(out_path,"/survival/"))) {
    dir.create(paste0(out_path,"/survival/"), recursive=TRUE)
}

setwd(paste0(out_path,"/survival/"))
## keep mapping for the aml samples
mapping_aml <- metadata_groups
mapping_aml$surv <- gsub("^0*", "",rownames(mapping_aml))
mapping_aml$surv <- gsub("[A-Z]?_.*", "", mapping_aml$surv)

## read in survival data
aml_surv <- read.csv("/data/metadata/Survival_AML_bulk_all.csv")
aml_surv_keep <- subset(aml_surv,aml_surv$ID %in% mapping_aml$surv)
survdata <- aml_surv_keep[,c("ID","OS.time","OS","DFI", "DFI.time")]

df.for.surv <- merge(survdata, mapping_aml, by.x="ID", by.y="surv")

## OS
df.for.surv$status.of.survival <- df.for.surv[,"OS"]
df.for.surv$OS.months <- df.for.surv[,"OS.time"]

## top 25th and bottom 25the percentile
df.for.surv$grp <- ifelse(df.for.surv$LSC_common_zscore >= quantile(df.for.surv$LSC_common_zscore ,0.75),"high",
                        ifelse(df.for.surv$LSC_common_zscore  <= quantile(df.for.surv$LSC_common_zscore ,0.25),"low",NA))
run_km(df.for.surv,paste0("Fig4c_", opname, "_quantile2575_OS"),"grp")

## LSC17 score survival in our cohort
df.for.surv$grp <- df.for.surv$LSC_high_low
run_km(df.for.surv,paste0("SupplFig6c_", opname, "_LSC17_high_low_OS"),"grp")


#DFI
df.for.surv$status.of.survival <- df.for.surv[,"DFI"]
df.for.surv$OS.months <- df.for.surv[,"DFI.time"]

## top 25th and bottom 25the percentile
df.for.surv$grp <- ifelse(df.for.surv$LSC_common_zscore >= quantile(df.for.surv$LSC_common_zscore ,0.75),"high",
                        ifelse(df.for.surv$LSC_common_zscore  <= quantile(df.for.surv$LSC_common_zscore ,0.25),"low",NA))
run_km(df.for.surv,paste0("Fig4d_", opname, "_quantile2575_DFI"),"grp")

# ## LSC17 score survival in our cohort
# df.for.surv$grp <- df.for.surv$LSC_high_low
# run_km(df.for.surv,paste0("SupplFig6e_", opname, "_LSC17_high_low_DFI"),"grp")



