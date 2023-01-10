# [CCBR] Filter Rows by Column Value (b47008ec-e06e-4540-bbff-bac339040af4): v12
Meta_VB013_Outlier_rem_noBASAL <- function(Meta_VB013_Outlier_rem) {
    library(tidyverse)
    tab <- Meta_VB013_Outlier_rem
    if(TRUE){
        x=paste0("^", c("BASAL"), "$", collapse = "|")    
    }
    else{
        x=paste(c("BASAL"),collapse="|")
    }
    if(FALSE){
    tab <- dplyr::filter(tab,grepl(x,TIMEPOINT,ignore.case = TRUE) )
    }
    else{
    tab <- dplyr::filter(tab,!grepl(x,TIMEPOINT,ignore.case = TRUE) )    
    }
    return(tab)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_Meta_VB013_Outlier_rem_noBASAL.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Meta_VB013_Outlier_rem<-readRDS(paste0(rds_output,"/var_Meta_VB013_Outlier_rem.rds"))
Input_is_Seurat_count <- 0
for(item in var_Meta_VB013_Outlier_rem){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Meta_VB013_Outlier_rem<-as.data.frame(var_Meta_VB013_Outlier_rem)}else{var_Meta_VB013_Outlier_rem <- var_Meta_VB013_Outlier_rem}
invisible(graphics.off())
var_Meta_VB013_Outlier_rem_noBASAL<-Meta_VB013_Outlier_rem_noBASAL(var_Meta_VB013_Outlier_rem)
invisible(graphics.off())
saveRDS(var_Meta_VB013_Outlier_rem_noBASAL, paste0(rds_output,"/var_Meta_VB013_Outlier_rem_noBASAL.rds"))
