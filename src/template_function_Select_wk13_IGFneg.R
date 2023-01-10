# [CCBR] Filter Rows by Column Value (b47008ec-e06e-4540-bbff-bac339040af4): v12
Select_wk13_IGFneg <- function(Select_IGFneg_VB013_P181) {
    library(tidyverse)
    tab <- Select_IGFneg_VB013_P181
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

print("template_function_Select_wk13_IGFneg.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Select_IGFneg_VB013_P181<-readRDS(paste0(rds_output,"/var_Select_IGFneg_VB013_P181.rds"))
Input_is_Seurat_count <- 0
for(item in var_Select_IGFneg_VB013_P181){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Select_IGFneg_VB013_P181<-as.data.frame(var_Select_IGFneg_VB013_P181)}else{var_Select_IGFneg_VB013_P181 <- var_Select_IGFneg_VB013_P181}
invisible(graphics.off())
var_Select_wk13_IGFneg<-Select_wk13_IGFneg(var_Select_IGFneg_VB013_P181)
invisible(graphics.off())
saveRDS(var_Select_wk13_IGFneg, paste0(rds_output,"/var_Select_wk13_IGFneg.rds"))
