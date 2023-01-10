Short_Gene_DEG_IGFneg_Age <- function(DEG_Age_IGF_neg) {
    library(limma)
    library(stringr)
    df <- dplyr::collect(DEG_Age_IGF_neg)
    df$Gene <- unlist(lapply(df$Gene,function(x) strsplit2(x,":")[[1]]))
    df$Gene <- unlist(lapply(df$Gene,function(x) strsplit2(x,"_")[[1]]))
    return(df)
}

print("template_function_Short_Gene_DEG_IGFneg_Age.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_DEG_Age_IGF_neg<-readRDS(paste0(rds_output,"/var_DEG_Age_IGF_neg.rds"))
Input_is_Seurat_count <- 0
for(item in var_DEG_Age_IGF_neg){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_DEG_Age_IGF_neg<-as.data.frame(var_DEG_Age_IGF_neg)}else{var_DEG_Age_IGF_neg <- var_DEG_Age_IGF_neg}
invisible(graphics.off())
var_Short_Gene_DEG_IGFneg_Age<-Short_Gene_DEG_IGFneg_Age(var_DEG_Age_IGF_neg)
invisible(graphics.off())
saveRDS(var_Short_Gene_DEG_IGFneg_Age, paste0(rds_output,"/var_Short_Gene_DEG_IGFneg_Age.rds"))
