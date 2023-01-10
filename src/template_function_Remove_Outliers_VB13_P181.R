# [CCBR] Filter Rows by Column Value (b47008ec-e06e-4540-bbff-bac339040af4): v12
Remove_Outliers_VB13_P181 <- function(Filtered_Study_VB013_P181) {
    library(tidyverse)
    tab <- Filtered_Study_VB013_P181
    if(TRUE){
        x=paste0("^", c("T112_wk13","H24G_wk13"), "$", collapse = "|")    
    }
    else{
        x=paste(c("T112_wk13","H24G_wk13"),collapse="|")
    }
    if(FALSE){
    tab <- dplyr::filter(tab,grepl(x,sampleID,ignore.case = TRUE) )
    }
    else{
    tab <- dplyr::filter(tab,!grepl(x,sampleID,ignore.case = TRUE) )    
    }
    return(tab)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_Remove_Outliers_VB13_P181.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Filtered_Study_VB013_P181<-readRDS(paste0(rds_output,"/var_Filtered_Study_VB013_P181.rds"))
Input_is_Seurat_count <- 0
for(item in var_Filtered_Study_VB013_P181){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Filtered_Study_VB013_P181<-as.data.frame(var_Filtered_Study_VB013_P181)}else{var_Filtered_Study_VB013_P181 <- var_Filtered_Study_VB013_P181}
invisible(graphics.off())
var_Remove_Outliers_VB13_P181<-Remove_Outliers_VB13_P181(var_Filtered_Study_VB013_P181)
invisible(graphics.off())
saveRDS(var_Remove_Outliers_VB13_P181, paste0(rds_output,"/var_Remove_Outliers_VB13_P181.rds"))
