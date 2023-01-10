# [CCBR] Filter Rows by Column Value (b47008ec-e06e-4540-bbff-bac339040af4): v12
Meta_P181_Outlier_rem <- function(Meta_P181) {
    library(tidyverse)
    tab <- Meta_P181
    if(TRUE){
        x=paste0("^", c("T112_wk13"), "$", collapse = "|")    
    }
    else{
        x=paste(c("T112_wk13"),collapse="|")
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

print("template_function_Meta_P181_Outlier_rem.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Meta_P181<-readRDS(paste0(rds_output,"/var_Meta_P181.rds"))
Input_is_Seurat_count <- 0
for(item in var_Meta_P181){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Meta_P181<-as.data.frame(var_Meta_P181)}else{var_Meta_P181 <- var_Meta_P181}
invisible(graphics.off())
var_Meta_P181_Outlier_rem<-Meta_P181_Outlier_rem(var_Meta_P181)
invisible(graphics.off())
saveRDS(var_Meta_P181_Outlier_rem, paste0(rds_output,"/var_Meta_P181_Outlier_rem.rds"))
