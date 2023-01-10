# [CCBR] Filter Rows by Column Value (b47008ec-e06e-4540-bbff-bac339040af4): v12
Select_IGFneg_VB013_P181 <- function(Remove_Outliers_VB13_P181) {
    library(tidyverse)
    tab <- Remove_Outliers_VB13_P181
    if(TRUE){
        x=paste0("^", c("0"), "$", collapse = "|")    
    }
    else{
        x=paste(c("0"),collapse="|")
    }
    if(TRUE){
    tab <- dplyr::filter(tab,grepl(x,IGF_1,ignore.case = TRUE) )
    }
    else{
    tab <- dplyr::filter(tab,!grepl(x,IGF_1,ignore.case = TRUE) )    
    }
    return(tab)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_Select_IGFneg_VB013_P181.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Remove_Outliers_VB13_P181<-readRDS(paste0(rds_output,"/var_Remove_Outliers_VB13_P181.rds"))
Input_is_Seurat_count <- 0
for(item in var_Remove_Outliers_VB13_P181){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Remove_Outliers_VB13_P181<-as.data.frame(var_Remove_Outliers_VB13_P181)}else{var_Remove_Outliers_VB13_P181 <- var_Remove_Outliers_VB13_P181}
invisible(graphics.off())
var_Select_IGFneg_VB013_P181<-Select_IGFneg_VB013_P181(var_Remove_Outliers_VB13_P181)
invisible(graphics.off())
saveRDS(var_Select_IGFneg_VB013_P181, paste0(rds_output,"/var_Select_IGFneg_VB013_P181.rds"))
