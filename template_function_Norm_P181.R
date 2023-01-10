Norm_P181 <- function(DEG_Age_IGF_neg, Meta_P181_Outlier_rem) {
    Meta <- Meta_P181_Outlier_rem
    DEG <- DEG_Age_IGF_neg
    
    library (tidyr)

    KeepColumns <- Meta$sampleID
    NormCounts <- DEG[c("Gene", KeepColumns)]

#    within(NormCounts, Gene<-data.frame(do.call('rbind', strsplit(as.character(Gene), '_', fixed=TRUE))))

   
    NormCounts <- separate(data = NormCounts, col = Gene, into = c("Gene", NA), sep = "_")

    return(NormCounts)
}

print("template_function_Norm_P181.R #########################################################################")
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
var_Meta_P181_Outlier_rem<-readRDS(paste0(rds_output,"/var_Meta_P181_Outlier_rem.rds"))
Input_is_Seurat_count <- 0
for(item in var_Meta_P181_Outlier_rem){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Meta_P181_Outlier_rem<-as.data.frame(var_Meta_P181_Outlier_rem)}else{var_Meta_P181_Outlier_rem <- var_Meta_P181_Outlier_rem}
invisible(graphics.off())
var_Norm_P181<-Norm_P181(var_DEG_Age_IGF_neg,var_Meta_P181_Outlier_rem)
invisible(graphics.off())
saveRDS(var_Norm_P181, paste0(rds_output,"/var_Norm_P181.rds"))
