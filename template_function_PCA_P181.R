# [CCBR] PCA Plot for RNA-Seq Data (1b7be3e3-4e2f-4390-b139-e114a699fac1): v126
PCA_P181 <- function(P181_LCF, Meta_P181_Outlier_rem){
    # image: png
    suppressMessages(library(ggplot2))
    suppressMessages(library(tidyverse))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(stringr))
    suppressMessages(library(RCurl))

    samples_to_include = c("R957_wk13","R958_wk13","R960_wk13","R962_wk13","R963_wk13","R966_wk13","R968_wk13","R976_wk13","T112_wk13","T114_wk13","T115_wk13","T121_wk13","T123_wk13","T125_wk13","T126_wk13","T127_wk13","T128_wk13")
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]

    P181_LCF %>% dplyr::select(append("Gene", samples_to_include)) -> spark.df
    
    spark.df %>% dplyr::collect() -> edf.orig

    Sampinfo <- dplyr::collect(Meta_P181_Outlier_rem)
    # cell <- Sampinfo$sampleID
    #Sampinfo = Sampinfo[, colSums(is.na(Sampinfo)) == 0]
    rownames(Sampinfo) <- Sampinfo$sampleID
    Sampinfo <- Sampinfo[match(colnames(edf.orig[,-1]), Sampinfo$sampleID), ]
    #Sampinfo = na.omit(Sampinfo)
    Sampinfo = Sampinfo[complete.cases(Sampinfo[, "sampleID"]),]
    print(paste0("Total number of genes included: ", nrow(edf.orig)))

    edf <- edf.orig[,match(Sampinfo$sampleID,colnames(edf.orig))]
    tedf <- t(edf)
    colnames(tedf) <- edf.orig[,1]
    tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
    tedf <- tedf[, apply(tedf, 2, var) != 0]
    pca <- prcomp(tedf, scale.=T)
    
    
    pca.df <- as.data.frame(pca$x) %>% dplyr::select(PC1, PC2)
    pca.df$celltype <- Sampinfo$ADJUVANT
    pca.df$sample <- Sampinfo$sampleID
    
    # Manual changes to sample names
    replacements = c("")

    plotcolors <- c("darkred","purple3","cadetblue","coral","deeppink","darkblue","darkgoldenrod","darkolivegreen3")
    if (length(unique(Sampinfo$Group)) > length(plotcolors)) {
        plotcolors <- c(plotcolors, rep("black", length(unique(Sampinfo$Group)) - length(plotcolors)))
    }

    if (!is.null(replacements)) {
        if (replacements != c("")) {
            for (x in replacements) {
                old <- strsplit(x, ": ?")[[1]][1]
                new <- strsplit(x, ": ?")[[1]][2]
                pca.df$sample <- ifelse(pca.df$sample==old, new, pca.df$sample)
            }
        }
    }
    perc.var <- (pca$sdev^2/sum(pca$sdev^2))*100
    perc.var <- formatC(perc.var,format = "g",digits=4)
    pc.x.lab <- paste0("PC1 ", perc.var[1],"%")
    pc.y.lab <- paste0("PC2 ", perc.var[2],"%")
    
    labelpos <- pca.df
    labelpos$mean_y <- pca.df$PC2+2
    labelpos$mean_x <- pca.df$PC1+2
    
    if (TRUE){
    g <- ggplot(pca.df, aes(x=PC1, y=PC2)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(legend.position="top") +
      geom_point(aes(color=pca.df$celltype), size=5) +
      geom_text(data=labelpos, aes(x=labelpos$mean_x, y=labelpos$mean_y, 
          label=sample, color=celltype, vjust="inward", hjust="inward"), size=5, show.legend=FALSE) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
      scale_colour_manual(values = plotcolors) +
      xlab(pc.x.lab) + ylab(pc.y.lab)
    }

    else{
    g <- ggplot(pca.df, aes(x=PC1, y=PC2)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(legend.position="top") +
      geom_point(aes(color=pca.df$celltype), size=5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
      scale_colour_manual(values = plotcolors) +
      xlab(pc.x.lab) + ylab(pc.y.lab)    
    }

    print(g)

    rownames(edf)=edf.orig[,1]
    edf.df = as.data.frame(edf)
    edf.df %>% rownames_to_column("Gene") -> edf.df
    return(edf.df)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_PCA_P181.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_P181_LCF<-readRDS(paste0(rds_output,"/var_P181_LCF.rds"))
Input_is_Seurat_count <- 0
for(item in var_P181_LCF){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_P181_LCF<-as.data.frame(var_P181_LCF)}else{var_P181_LCF <- var_P181_LCF}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Meta_P181_Outlier_rem<-readRDS(paste0(rds_output,"/var_Meta_P181_Outlier_rem.rds"))
Input_is_Seurat_count <- 0
for(item in var_Meta_P181_Outlier_rem){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Meta_P181_Outlier_rem<-as.data.frame(var_Meta_P181_Outlier_rem)}else{var_Meta_P181_Outlier_rem <- var_Meta_P181_Outlier_rem}
invisible(graphics.off())
var_PCA_P181<-PCA_P181(var_P181_LCF,var_Meta_P181_Outlier_rem)
invisible(graphics.off())
saveRDS(var_PCA_P181, paste0(rds_output,"/var_PCA_P181.rds"))
