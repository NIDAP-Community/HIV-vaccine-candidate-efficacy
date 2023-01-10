# [CCBR] Filter Low Count Genes For Bulk RNAseq Data (59e7ee62-1715-4276-b370-e8395168f9d8): v57
Filter_low_count_IGFneg <- function(Normalized_Data, Select_wk13_IGFneg) {
    suppressMessages(library(limma))
    suppressMessages(library(tidyverse))
    suppressMessages(library(edgeR))
    suppressMessages(library(ggplot2))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(stringr))
    suppressMessages(library(RCurl))
    suppressMessages(library(reshape))
    suppressMessages(library(gridExtra))

    spark.df = Normalized_Data
    samples_to_include = c("H02V_BASAL","H02V_wk13","H21Y_wk13","H22D_wk13","H22E_BASAL","H22E_wk13","H22F_wk13","H22H_BASAL","H24C_BASAL","H24D_wk13","H24G_wk13","H24H_wk13","H28T_wk13","H28V_wk13","H28Y_wk13","H29N_BASAL","H30G_BASAL","H34F_BASAL","H35H_wk13","H35R_wk13","R957_wk13","R958_wk13","R960_wk13","R962_wk13","R963_wk13","R966_wk13","R968_wk13","R976_wk13","T112_wk13","T114_wk13","T115_wk13","T121_wk13","T123_wk13","T125_wk13","T126_wk13","T127_wk13","T128_wk13")
    samples_to_include <- samples_to_include[samples_to_include != "mir_id"]
    spark.df %>% dplyr::select(append("mir_id", samples_to_include)) -> spark.df

    df <- drop_na(spark.df)
    df <- dplyr::collect(df)
    df %>% dplyr::group_by(mir_id) %>% summarize_all(sum) -> df
    print(paste0("Number of genes before filtering: ", nrow(df)))
    df$isexpr1 <- rowSums(cpm(as.matrix(df[, -1])) > 1) >= 1
    df <- as.data.frame(df[df$isexpr1, ])
    df$isexpr1 <- NULL
    
    if (TRUE) {
        rownames(df) <- df$mir_id
        df$mir_id <- NULL
        counts <- cpm(as.matrix(df)) > 1 # boolean matrix
        tcounts <- as.data.frame(t(counts))
        colnum <- dim(counts)[1] # number of genes
        metadata <- dplyr::collect(Select_wk13_IGFneg)
        rownames(metadata) <- metadata$sampleID
        tcounts <- merge(metadata["STUDY"], tcounts, by="row.names")
        tcounts$Row.names <- NULL
        melted <- melt(tcounts, id.vars="STUDY")
        tcounts.tot <- summarise(group_by(melted, STUDY, variable), sum=sum(value))
        tcounts.tot %>% tidyr::spread(variable, sum) -> tcounts.group
        colSums(tcounts.group[(1:colnum+1)]>=3) >= 1 -> tcounts.keep
        df <- df[tcounts.keep, ]
        df %>% rownames_to_column("mir_id") -> df
    }

    colnames(df)[colnames(df)=="mir_id"] <- "Gene"

    print(paste0("Number of genes after filtering: ", nrow(df)))
    
    df -> edf.orig

    Sampinfo <- dplyr::collect(Select_wk13_IGFneg)
    # cell <- Sampinfo$sampleID
    Sampinfo = Sampinfo[, colSums(is.na(Sampinfo)) == 0]
    rownames(Sampinfo) <- Sampinfo$sampleID
    Sampinfo <- Sampinfo[match(colnames(edf.orig[,-1]), Sampinfo$sampleID), ]
    Sampinfo = Sampinfo[!is.na(Sampinfo["sampleID"]),]
    print(paste0("Total number of genes included: ", nrow(edf.orig)))

    edf <- edf.orig[,match(Sampinfo$sampleID,colnames(edf.orig))]
    tedf <- t(edf)
    colnames(tedf) <- edf.orig[,1]
    tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
    tedf <- tedf[, apply(tedf, 2, var) != 0]
    pca <- prcomp(tedf, scale.=T)
    
    pca.df <- as.data.frame(pca$x) %>% dplyr::select(PC1, PC2)
    pca.df$celltype <- Sampinfo$STUDY
    pca.df$sample <- Sampinfo$sampleID
    
    # Manual changes to sample names
    replacements = c()

    plotcolors <- c("aquamarine3","salmon1","lightskyblue3","plum3","darkolivegreen3","goldenrod1","burlywood2","gray70","firebrick2","steelblue","palegreen4","orchid4","darkorange1","yellow","sienna","palevioletred1","gray60","cyan4","darkorange3","mediumpurple3","violetred2","olivedrab","darkgoldenrod2","darkgoldenrod","gray40","palegreen3","thistle3","khaki1","deeppink2","chocolate3","paleturquoise3","wheat1","lightsteelblue","salmon","sandybrown","darkolivegreen2","thistle2","gray85","orchid3","darkseagreen1","lightgoldenrod1","lightskyblue2","dodgerblue3","darkseagreen3","forestgreen","lightpink2","mediumpurple4","lightpink1","thistle","navajowhite","lemonchiffon","bisque2","mistyrose","gray95","lightcyan3","peachpuff2","lightsteelblue2","lightyellow2","moccasin","antiquewhite2","gray80","lightgrey")
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
      geom_point(aes(color=pca.df$celltype), size=1) +
      geom_text(data=labelpos, aes(x=labelpos$mean_x, y=labelpos$mean_y, 
          label=sample, color=celltype, vjust="inward", hjust="inward"), size=3, show.legend=FALSE) +
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
      geom_point(aes(color=pca.df$celltype), size=1) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
      scale_colour_manual(values = plotcolors) +
      xlab(pc.x.lab) + ylab(pc.y.lab)    
    }
    par(mfrow = c(2,1))
   
    # samples_to_include = c("H02V_BASAL","H02V_wk13","H21Y_wk13","H22D_wk13","H22E_BASAL","H22E_wk13","H22F_wk13","H22H_BASAL","H24C_BASAL","H24D_wk13","H24G_wk13","H24H_wk13","H28T_wk13","H28V_wk13","H28Y_wk13","H29N_BASAL","H30G_BASAL","H34F_BASAL","H35H_wk13","H35R_wk13","R957_wk13","R958_wk13","R960_wk13","R962_wk13","R963_wk13","R966_wk13","R968_wk13","R976_wk13","T112_wk13","T114_wk13","T115_wk13","T121_wk13","T123_wk13","T125_wk13","T126_wk13","T127_wk13","T128_wk13")
    # samples_to_include <- samples_to_include[samples_to_include != "mir_id"]
    df.filt <- df %>% select(samples_to_include)
    Sampinfo <- dplyr::collect(Select_wk13_IGFneg)
        # cell <- Sampinfo$SampleName
        rownames(Sampinfo) <- Sampinfo$sampleID
        Sampinfo <- Sampinfo[match(colnames(df.filt[,-1]), Sampinfo$sampleID), ]
        Sampinfo =  Sampinfo[!is.na(Sampinfo["sampleID"]),]
        print(paste0("Total number of samples included: ", nrow(Sampinfo)))
    df.filt <- df.filt[,match(Sampinfo$sampleID,colnames(df.filt))]
    rownames(df.filt)=df[,1]
    df.filt %>% rownames_to_column("mir_id") -> df.filt
    df.m <- melt(df.filt,id.vars=c("mir_id"))
    if(TRUE){
    df.m %>% mutate(value = log2(value+0.5)) -> df.m
    }
    n <- 40
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    qual_col_pals = qual_col_pals[c(7,6,2,1,8,3,4,5),]
    cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

    if(FALSE){
        xmin = -1
        xmax = 1
    }
    else{
        xmin = min(df.m$value)
        xmax = max(df.m$value)
    }

    if(FALSE){
    df.m %>% mutate(colgroup = Sampinfo[variable,"STUDY"]) -> df.m
    g2=ggplot(df.m, aes(x=value, group=variable)) + 
        geom_density(aes(colour = colgroup, linetype = colgroup))+
        xlab(NULL) + ylab(NULL) +
         theme_bw() +
        theme(legend.position='top', legend.text = element_text(size = 5),legend.title=element_blank()) + #scale_x_log10() + 
        xlim(xmin,xmax) +
        scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),40)) + 
        scale_colour_manual(values=cols) 
    }
    else{
        df.m$variable = Sampinfo[df.m$variable,"sampleID"]
     n=length(unique(df.m$variable))
     if(n>160){
     m=ceiling(n/4)
     n=m*4
     color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
     cols=sample(color, n)}
     else{
         m=40
     }
        g2=ggplot(df.m, aes(x=value, group=variable)) + 
        geom_density(aes(colour = variable, linetype = variable) ) + 
        xlab(NULL) + ylab(NULL) +
        theme_bw() +
        theme(legend.position='top', legend.title=element_blank(), legend.text = element_text(size = 5)) + 
        xlim(xmin,xmax) +
        scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),m)) + 
        scale_colour_manual(values=cols)
    }
    require(gridExtra)
    plots = grid.arrange(g,g2, ncol=2, respect=TRUE)
    print(plots)

    return(df)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

print("template_function_Filter_low_count_IGFneg.R #########################################################################")
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
var_Select_wk13_IGFneg<-readRDS(paste0(rds_output,"/var_Select_wk13_IGFneg.rds"))
Input_is_Seurat_count <- 0
for(item in var_Select_wk13_IGFneg){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Select_wk13_IGFneg<-as.data.frame(var_Select_wk13_IGFneg)}else{var_Select_wk13_IGFneg <- var_Select_wk13_IGFneg}
invisible(graphics.off())
var_Filter_low_count_IGFneg<-Filter_low_count_IGFneg(var_Normalized_Data,var_Select_wk13_IGFneg)
invisible(graphics.off())
saveRDS(var_Filter_low_count_IGFneg, paste0(rds_output,"/var_Filter_low_count_IGFneg.rds"))
