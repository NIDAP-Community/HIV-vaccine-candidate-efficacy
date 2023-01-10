# [CCBR] Similarity Heatmap with ggplot (1bc3e623-b220-4ef9-9984-7a61b774be6f): v1
Similarity_P181 <- function(P181_LCF,MacacaMetadata) {
suppressMessages(library(amap))
suppressMessages(library(lattice))
suppressMessages(library(gplots))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(reshape2))
suppressMessages(library(ggdendro))
suppressMessages(library(scales))
suppressMessages(library(gtable))
df <- dplyr::collect(P181_LCF)
samples_to_include = c("R957_wk13","R958_wk13","R960_wk13","R962_wk13","R963_wk13","R966_wk13","R968_wk13","R976_wk13","T112_wk13","T114_wk13","T115_wk13","T121_wk13","T123_wk13","T125_wk13","T126_wk13","T127_wk13","T128_wk13")
samples_to_include <- samples_to_include[samples_to_include != "Gene"]
df <- df[,samples_to_include]

targetfile <- dplyr::collect(MacacaMetadata)
targetfile %>% filter(sampleID %in% samples_to_include) -> targetfile
rownames(targetfile) = targetfile$sampleID
df <- df[,match(targetfile$sampleID,colnames(df))]

  mat <- as.matrix(t(df)) 
  d=Dist(mat,method="correlation",diag=TRUE)
  m=as.matrix(d)
  h = hclust(d)
  dim(m)
  
  dd.col <- as.dendrogram(h)
  col.ord <- order.dendrogram(dd.col)
  
  dd.row <- as.dendrogram(h)
  row.ord <- order.dendrogram(dd.row)

  xx <- as.data.frame(m)[col.ord, row.ord]
  df <- as.data.frame(xx)
  colnames(df) <- attr(xx,c("names"))

  df$y.variable <- attr(xx,c("names"))
  df$y.variable <- with(df, factor(y.variable, 
        levels=y.variable, ordered=TRUE))

  mdf <- melt(df, id.vars="y.variable")
  
  ### Use ggdendro to extract dendrogram data ###
  ddata_x <- dendro_data(dd.row)
  ddata_y <- dendro_data(dd.col)
  
  ##Set up default theme:

    theme_none <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
)

  ### Create plot components ###    
  # Heatmap
  p1 <- ggplot(mdf,aes(x=variable, y=y.variable,fill=value)) + 
    geom_tile()+
    theme_none +
    scale_fill_gradient2(low = muted("blue4"), mid = "white", high = "red",
                          midpoint = median(mdf$value)) +
    #scale_fill_gradient(low="green", high="red") +
    scale_x_discrete(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0))
  
  #Add colorbars:
  vardf <- data.frame(var=targetfile[unique(mdf$variable),"censor"],
      xpos=1,ypos=c(1:length(unique(mdf$variable))))
  
  #Top
  p1a <- ggplot(vardf,aes(xpos,ypos)) + 
    geom_tile(aes(fill=as.character(var)),width=1) + 
    coord_flip() + theme_none + 
    scale_fill_brewer(type="qual",palette = "Set2") +
    #scale_fill_brewer() +
    scale_x_discrete(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0))

  #Colorbar (left):
  p1l <- ggplot(vardf,aes(xpos,ypos)) + 
    geom_tile(aes(fill=as.character(var)),width=1) + 
    #coord_flip() + 
    theme_none + 
    #theme(legend.position = "right", legend.title = element_blank()) +
    scale_fill_brewer(type="qual",palette = "Set2") +
    #scale_fill_brewer() +
    scale_x_discrete(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0))
  p1l
  gp1l <- ggplotGrob(p1l)
  
  p1aleg <- ggplot(vardf,aes(xpos,ypos)) + 
    geom_tile(aes(fill=as.character(var)),width=1) + 
    coord_flip() + theme_none + 
    theme(legend.position = "right", legend.text = element_text(size=20), legend.title = element_blank()) +
    scale_fill_brewer(type="qual",palette = "Set2") +
    #scale_fill_brewer() +
    scale_x_discrete(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0))
  gp1aleg <- ggplotGrob(p1aleg)
  l <- gtable_filter(gp1aleg, 'guide-box', trim=F)

  
  p1b <- ggplot(mdf,aes(x=variable, y=y.variable)) +
    theme_none +
    theme(axis.text.x = element_text(size=12,angle = 90, vjust = 0.5, hjust=1)) +
    scale_x_discrete(expand = expand_scale(add = .6)) + 
    scale_y_discrete(expand = c(0, 0))
  gp1b <- ggplotGrob(p1b)
  s <- gtable_filter(gp1b, 'axis-b|xlab', trim=F)  # use trim depending on need
  #grid.newpage()
  #grid.draw(s)
  
   #Add y-axis:
   p1c <- ggplot(mdf,aes(x=variable, y=y.variable)) +
    theme_none +
    theme(axis.text.y = element_text(size=8, vjust = 0.5, hjust=1)) +
    scale_x_discrete(expand = c(0,0)) + 
    scale_y_discrete(expand = c(0, 0))
  gp1c <- ggplotGrob(p1c)
  t <- gtable_filter(gp1c, 'axis-l|ylab', trim=F) 

  # Dendrogram 1 (top)
  p2 <- ggplot(segment(ddata_x)) + 
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
    theme_none +
    scale_x_continuous(expand = expand_scale(add = .6)) +
    scale_y_continuous(expand = c(0, 0))

  gp1 <- ggplotGrob(p1)
  gp2 <- ggplotGrob(p2)
  gp1a <- ggplotGrob(p1a)

  lay <- rbind(c(NA,NA,1,NA),
             c(NA,NA,2,NA),
             c(3,4,5,6),
             c(NA,NA,7,NA))
  gcomb = grid.arrange(gp2,gp1a,t,gp1l,gp1,l,s, ncol=4, nrow=4, heights=c(0.5,0.25,3,0.5),widths=c(0.25,0.25,4,0.25),layout_matrix=lay)
  #gcomb = grid.arrange(gp2,gp1a,gp1,s,ncol=1,nrow=4,heights=c(0.5,0.25,3,0.5))
  #gcomb = arrangeGrob(gp2,gp1a,gp1,s,ncol=1,nrow=4, heights=c(0.5,0.25,3,0.5))
  #gcomb[["heights"]] <- c(0.5,0.25,3,0.5)
  #print(gcomb[["heights"]])

  #gcomb2 = grid.arrange(gcomb,l,ncol=2,nrow=1,widths=c(4,0.25))
  print(gcomb)

return(P181_LCF)   
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

print("template_function_Similarity_P181.R #########################################################################")
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
var_MacacaMetadata<-readRDS(paste0(rds_output,"/var_MacacaMetadata.rds"))
Input_is_Seurat_count <- 0
for(item in var_MacacaMetadata){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_MacacaMetadata<-as.data.frame(var_MacacaMetadata)}else{var_MacacaMetadata <- var_MacacaMetadata}
invisible(graphics.off())
var_Similarity_P181<-Similarity_P181(var_P181_LCF,var_MacacaMetadata)
invisible(graphics.off())
saveRDS(var_Similarity_P181, paste0(rds_output,"/var_Similarity_P181.rds"))
