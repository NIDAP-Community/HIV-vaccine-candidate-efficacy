Norm_VB013 <- function(DEG_Age_IGF_neg, Meta_VB013_Outlier_rem_noBASAL) {
    Meta <- Meta_VB013_Outlier_rem_noBASAL
    DEG <- DEG_Age_IGF_neg
    
    library (tidyr)

    KeepColumns <- Meta$sampleID
    NormCounts <- DEG[c("Gene", KeepColumns)]

#    within(NormCounts, Gene<-data.frame(do.call('rbind', strsplit(as.character(Gene), '_', fixed=TRUE))))

   
    NormCounts <- separate(data = NormCounts, col = Gene, into = c("Gene", NA), sep = "_")

    return(NormCounts)
}

print("template_function_Norm_VB013.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_DEG_Age_IGF_neg<-readRDS(paste0(rds_output,"/var_DEG_Age_IGF_neg.rds"))
Input_is_Seurat_count <- 0
for(item in var_DEG_Age_IGF_neg){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_DEG_Age_IGF_neg<-as.data.frame(var_DEG_Age_IGF_neg)}else{var_DEG_Age_IGF_neg <- var_DEG_Age_IGF_neg}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Meta_VB013_Outlier_rem_noBASAL<-readRDS(paste0(rds_output,"/var_Meta_VB013_Outlier_rem_noBASAL.rds"))
Input_is_Seurat_count <- 0
for(item in var_Meta_VB013_Outlier_rem_noBASAL){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Meta_VB013_Outlier_rem_noBASAL<-as.data.frame(var_Meta_VB013_Outlier_rem_noBASAL)}else{var_Meta_VB013_Outlier_rem_noBASAL <- var_Meta_VB013_Outlier_rem_noBASAL}
invisible(graphics.off())
var_Norm_VB013<-Norm_VB013(var_DEG_Age_IGF_neg,var_Meta_VB013_Outlier_rem_noBASAL)
invisible(graphics.off())
saveRDS(var_Norm_VB013, paste0(rds_output,"/var_Norm_VB013.rds"))
