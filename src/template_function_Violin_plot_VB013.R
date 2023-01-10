Violin_plot_VB013 <- function(Normalized_Data, MacacaMetadata) {

    suppressMessages(library(tidyr))
    suppressMessages(library(reshape))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(ggplot2))

    df <- dplyr::collect(Normalized_Data)
    samples <- c("mir_id","H02V_BASAL","H02V_wk13","H21Y_wk13","H22D_wk13","H22E_BASAL","H22E_wk13","H22F_wk13","H22H_BASAL","H24C_BASAL","H24D_wk13","H24G_wk13","H24H_wk13","H28T_wk13","H28V_wk13","H28Y_wk13","H29N_BASAL","H30G_BASAL","H34F_BASAL","H35H_wk13","H35R_wk13","R957_wk13","R958_wk13","R960_wk13","R962_wk13","R963_wk13","R966_wk13","R968_wk13","R976_wk13","T112_wk13","T114_wk13","T115_wk13","T121_wk13","T123_wk13","T125_wk13","T126_wk13","T127_wk13","T128_wk13")
    df <- df[,samples] 
    
    targetfile <- dplyr::collect(MacacaMetadata)
    rownames(targetfile) = targetfile$sampleID
    targetfile <- targetfile %>% filter(sampleID %in% samples) %>% arrange(TIMEPOINT)
    df = df[,match(targetfile$sampleID,colnames(df))]
    df.m = melt(df)

    
    if(TRUE){
    df.m$value = log2(df.m$value + 0.5)
    }

    df.m$group = targetfile[df.m$variable,"TIMEPOINT"]
    min.num = min(df.m$value)-(1.1*abs(min(df.m$value)))
    max.num = max(df.m$value)*1.1
    #min.num = 0
    #max.num = 5
    col1=brewer.pal(8, "Set3")[-2]
    
    g = ggplot(df.m, aes(x=variable, y=value)) +
    #ggtitle(gene) +
    labs(y = "Expression Values",size=2) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 20, face="italic"),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14), 
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_violin(aes(fill=as.factor(group))) + 
    scale_fill_manual(values = col1) +
    coord_cartesian(ylim=c(min.num,max.num)) +
    #labs(colour = n, y=m) +
    geom_jitter(height = 0, width = 0.1, size = 0.1) +
    ylim(0,10)

  g + geom_boxplot(width=0.05) #+
  #    coord_cartesian(ylim=c(min.num,max.num))

print(g)
return(df)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_Violin_plot_VB013.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Normalized_Data<-readRDS(paste0(rds_output,"/var_Normalized_Data.rds"))
Input_is_Seurat_count <- 0
for(item in var_Normalized_Data){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Normalized_Data<-as.data.frame(var_Normalized_Data)}else{var_Normalized_Data <- var_Normalized_Data}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_MacacaMetadata<-readRDS(paste0(rds_output,"/var_MacacaMetadata.rds"))
Input_is_Seurat_count <- 0
for(item in var_MacacaMetadata){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_MacacaMetadata<-as.data.frame(var_MacacaMetadata)}else{var_MacacaMetadata <- var_MacacaMetadata}
invisible(graphics.off())
var_Violin_plot_VB013<-Violin_plot_VB013(var_Normalized_Data,var_MacacaMetadata)
invisible(graphics.off())
saveRDS(var_Violin_plot_VB013, paste0(rds_output,"/var_Violin_plot_VB013.rds"))
