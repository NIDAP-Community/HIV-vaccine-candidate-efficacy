TRUNC_DEG_VB103 <- function(DEG_VB013) {
    library(limma)
    library(stringr)
    df <- dplyr::collect(DEG_VB013)
    df$Gene <- unlist(lapply(df$Gene,function(x) strsplit2(x,":")[[1]]))
    df$Gene <- unlist(lapply(df$Gene,function(x) strsplit2(x,"_")[[1]]))
    return(df)
}

print("template_function_TRUNC_DEG_VB103.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_DEG_VB013<-readRDS(paste0(rds_output,"/var_DEG_VB013.rds"))
Input_is_Seurat_count <- 0
for(item in var_DEG_VB013){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_DEG_VB013<-as.data.frame(var_DEG_VB013)}else{var_DEG_VB013 <- var_DEG_VB013}
invisible(graphics.off())
var_TRUNC_DEG_VB103<-TRUNC_DEG_VB103(var_DEG_VB013)
invisible(graphics.off())
saveRDS(var_TRUNC_DEG_VB103, paste0(rds_output,"/var_TRUNC_DEG_VB103.rds"))
