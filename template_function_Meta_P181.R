# [CCBR] Filter Rows by Column Value (b47008ec-e06e-4540-bbff-bac339040af4): v12
Meta_P181 <- function(MacacaMetadata) {
    library(tidyverse)
    tab <- MacacaMetadata
    if(TRUE){
        x=paste0("^", c("P181"), "$", collapse = "|")    
    }
    else{
        x=paste(c("P181"),collapse="|")
    }
    if(TRUE){
    tab <- dplyr::filter(tab,grepl(x,STUDY,ignore.case = TRUE) )
    }
    else{
    tab <- dplyr::filter(tab,!grepl(x,STUDY,ignore.case = TRUE) )    
    }
    return(tab)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_Meta_P181.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_MacacaMetadata<-readRDS(paste0(rds_output,"/var_MacacaMetadata.rds"))
Input_is_Seurat_count <- 0
for(item in var_MacacaMetadata){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_MacacaMetadata<-as.data.frame(var_MacacaMetadata)}else{var_MacacaMetadata <- var_MacacaMetadata}
invisible(graphics.off())
var_Meta_P181<-Meta_P181(var_MacacaMetadata)
invisible(graphics.off())
saveRDS(var_Meta_P181, paste0(rds_output,"/var_Meta_P181.rds"))
