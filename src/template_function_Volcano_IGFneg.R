# [CCBR] Volcano Plot (eeced39d-ed52-4b16-9847-4282971775e6): v62
Volcano_IGFneg <- function(Short_Gene_DEG_IGFneg_Age) {
    # image: png
    suppressMessages(library(ggplot2))
    suppressMessages(library(dplyr))
    suppressMessages(library(ggrepel))
    genesmat <- dplyr::collect(Short_Gene_DEG_IGFneg_Age)
    value_to_sort_the_output_dataset = "p-value"
    if (value_to_sort_the_output_dataset=="fold-change") {
        genesmat %>% arrange(desc(abs(`P181-VB013_logFC`))) -> genesmat
    } else if (value_to_sort_the_output_dataset=="p-value") {
        genesmat %>% arrange(`P181-VB013_pval`) -> genesmat
    }
    print(paste0("Total number of genes included in volcano plot: ", nrow(genesmat)))
    if (TRUE){
    ymax = ceiling(max(-log10(genesmat$`P181-VB013_pval`)))}
    else (ymax = 10)
    if (TRUE){
    xmax = ceiling(max(genesmat$`P181-VB013_logFC`))}
    else (xmax = 5)
    
    p <- ggplot(genesmat,
        aes(x = `P181-VB013_logFC`, y = -log10(`P181-VB013_pval`))) +
        theme_classic() +
        geom_point(
            color='black',
            size = 1) +
        geom_vline(xintercept=c(-1,1), color='red', alpha=1.0) + 
        geom_hline(yintercept=-log10(0.001), color='blue', alpha=1.0) +  
        geom_point(
            data = subset(genesmat, (`P181-VB013_pval` < 0.001)),
            color = 'lightgoldenrod2',
            size = 1) +
        geom_point(
            data = subset(genesmat, (`P181-VB013_pval` < 0.001 & abs(`P181-VB013_logFC`)>1)),
            color = 'red',
            size = 1) +
        geom_text_repel(
            data = genesmat[1:30,],
            label = genesmat$Gene[1:30],
            nudge_x = 0.2,
            nudge_y = 0.2,
            size = 4) +
        xlim(-xmax,xmax) +
        ylim(0,ymax) 
    print(p)

    return(genesmat)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_Volcano_IGFneg.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Short_Gene_DEG_IGFneg_Age<-readRDS(paste0(rds_output,"/var_Short_Gene_DEG_IGFneg_Age.rds"))
Input_is_Seurat_count <- 0
for(item in var_Short_Gene_DEG_IGFneg_Age){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Short_Gene_DEG_IGFneg_Age<-as.data.frame(var_Short_Gene_DEG_IGFneg_Age)}else{var_Short_Gene_DEG_IGFneg_Age <- var_Short_Gene_DEG_IGFneg_Age}
invisible(graphics.off())
var_Volcano_IGFneg<-Volcano_IGFneg(var_Short_Gene_DEG_IGFneg_Age)
invisible(graphics.off())
saveRDS(var_Volcano_IGFneg, paste0(rds_output,"/var_Volcano_IGFneg.rds"))
