# [CCBR] Voom Normalization and DEG Analysis (de953b9c-b7f3-49ac-b883-55070d759e5d): v149
DEG_VB013 <- function(VB013_LCF, Meta_VB013_Outlier_rem) {
    suppressMessages(library(limma))
    suppressMessages(library(tidyverse))
    suppressMessages(library(edgeR))
    suppressMessages(library(stringr))
    df <- VB013_LCF 

    if(make.names(colnames(df))!=colnames(df)){
        print("Error: The following counts matrix column names are not valid:\n")
        print(colnames(df)[make.names(colnames(df))!=colnames(df)])

        print("Likely causes are columns starting with numbers or other special characters eg spaces.")
        stop("Bad column names.")
    }

    samples_for_deg_analysis = c("H02V_BASAL","H02V_wk13","H21Y_wk13","H22D_wk13","H22E_BASAL","H22E_wk13","H22F_wk13","H22H_BASAL","H24C_BASAL","H24D_wk13","H24G_wk13","H24H_wk13","H28T_wk13","H28V_wk13","H28Y_wk13","H29N_BASAL","H30G_BASAL","H34F_BASAL","H35H_wk13","H35R_wk13")
    targetfile <- Meta_VB013_Outlier_rem
    contrasts_param<-c("wk13-BASAL")
    sample_name_column_param<-"sampleID"
    ordered_covariates=c("TIMEPOINT","animal")
    return_mean_and_sd_param<-TRUE
    contrast_variable_param<-"TIMEPOINT"
    opt_normalization_method<-"quantile"

    samples_for_deg_analysis <- samples_for_deg_analysis[samples_for_deg_analysis != "Gene"]
    samples_for_deg_analysis <- samples_for_deg_analysis[samples_for_deg_analysis != "GeneName"]

    df.m <- df[,samples_for_deg_analysis]

    gene_names <- NULL
    gene_names$GeneID <- df[,1]
    
    
    targetfile <- targetfile[match(colnames(df.m),targetfile[,sample_name_column_param]),]
    targetfile <- targetfile[rowSums(is.na(targetfile)) != ncol(targetfile), ]
    df.m <- df.m[,match(targetfile[,sample_name_column_param],colnames(df.m))]
    if(FALSE){
    x <- DGEList(counts=2^df.m, genes=gene_names)
    }
    else{
    x <- DGEList(counts=df.m, genes=gene_names)     
    }
    
    
    ordered_covariates=ordered_covariates[order(ordered_covariates!=contrast_variable_param)]

    dm.formula <- as.formula(paste("~0 +", paste(ordered_covariates, sep="+", collapse="+")))
    design=model.matrix(dm.formula, targetfile)
    #print(dm.formula)
    #print(design)
    colnames(design) <- str_replace_all(colnames(design), contrast_variable_param, "")
    
    if (opt_normalization_method %in% c("TMM","TMMwzp","RLE","upperquartile")){
    x <- calcNormFactors(x, method = opt_normalization_method) 
    rownames(x) <- x$genes$GeneID
    v <- voom_v2(x,design=design,normalize="none")
    }
    else {
    v <- voom_v2(x,design=design,normalize=opt_normalization_method,plot=TRUE)
    }
    
    rownames(v$E) <- v$genes$GeneID
    as.data.frame(v$E) %>% rownames_to_column("Gene") -> df.voom
    fit <- lmFit(v, design)
   
    cm <- makeContrasts(contrasts = contrasts_param, levels=design)
    
    fit2 <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit2)

    logFC = fit2$coefficients
    colnames(logFC)=paste(colnames(logFC),"logFC",sep="_")
    tstat = fit2$t
    colnames(tstat)=paste(colnames(tstat),"tstat",sep="_")
    FC = 2^fit2$coefficients
    FC = ifelse(FC<1,-1/FC,FC)
    colnames(FC)=paste(colnames(FC),"FC",sep="_")
    pvalall=fit2$p.value
    colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")
    pvaladjall=apply(pvalall,2,function(x) p.adjust(x,"BH"))
    colnames(pvaladjall)=paste(colnames(fit2$coefficients),"adjpval",sep="_")
    
    
    if(return_mean_and_sd_param ){
    tve <- t(v$E)
    mean.df <- as.data.frame(tve) %>% rownames_to_column("Sample") %>% mutate(group=targetfile[targetfile[,sample_name_column_param]==Sample,contrast_variable_param]) %>% group_by(group) %>% summarise_all(funs(mean)) %>% as.data.frame()
    mean.df[,-c(1,2)] %>% as.matrix() %>% t() -> mean
    colnames(mean) <- mean.df[,1]
    colnames(mean)=paste("Mean",colnames(mean),sep="_")
    colnames(mean) = gsub("\\.", "_", colnames(mean))
    sd.df <- as.data.frame(tve) %>% rownames_to_column("Sample") %>% mutate(group=targetfile[targetfile[,sample_name_column_param]==Sample,contrast_variable_param]) %>% group_by(group) %>% summarise_all(funs(sd)) %>% as.data.frame()
    sd.df[,-c(1,2)] %>% as.matrix() %>% t() -> sd
    colnames(sd) <- sd.df[,1]
    colnames(sd)=paste("SD",colnames(sd),sep="_")
    colnames(sd) = gsub("\\.", "_", colnames(sd))
   finalres=as.data.frame(cbind(mean,sd,v$E,FC, logFC, tstat, pvalall, pvaladjall)) 
}
   else{
    finalres=as.data.frame(cbind(v$E,FC, logFC, tstat, pvalall, pvaladjall))
    }
    finalres %>% rownames_to_column("Gene") -> finalres
    print(paste0("Total number of genes included: ", nrow(finalres)))
    
    call_me_alias<-colnames(finalres)
    colnames(finalres)<-gsub("\\(|\\)","",call_me_alias)
    sparkdffr<-(finalres)
   # print("Returning sparkR dataframe")
    return(sparkdffr) 
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

#voom function was modified to incrase dot size,line-width, x and y Label size, title label size, margin between tick, tick-values, and labels.

voom_v2 <- function (counts, design = NULL, lib.size = NULL, normalize.method = "none", 
                     span = 0.5, plot = FALSE, save.plot = FALSE, ...) 
{
  out <- list()
  if (is(counts, "DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) > 
        0) 
      design <- model.matrix(~group, data = counts$samples)
    if (is.null(lib.size)) 
      lib.size <- counts$samples$lib.size * counts$samples$norm.factors
    counts <- counts$counts
  }
  else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts, 
                                                         "ExpressionSet"))
    if (isExpressionSet) {
      if (length(Biobase::fData(counts))) 
        out$genes <- Biobase::fData(counts)
      if (length(Biobase::pData(counts))) 
        out$targets <- Biobase::pData(counts)
      counts <- Biobase::exprs(counts)
    }
    else {
      counts <- as.matrix(counts)
    }
  }
  n <- nrow(counts)
  if (n < 2L) 
    stop("Need at least two genes to fit a mean-variance trend")
  if (is.null(design)) {
    design <- matrix(1, ncol(counts), 1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  if (is.null(lib.size)) 
    lib.size <- colSums(counts)
  y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  y <- normalizeBetweenArrays(y, method = normalize.method)
  fit <- lmFit(y, design, ...)
  if (is.null(fit$Amean)) 
    fit$Amean <- rowMeans(y, na.rm = TRUE)
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  if (plot) {
    par(mgp=c(7,2.5,0),mar=c(9,11,6,4)+.1)
    #The 'mar' argument of 'par' sets the width of the margins in the order: 'bottom', 'left', 'top', 'right'. 
    plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", 
         pch = 16, cex = 1.25, cex.lab=4.5,cex.axis=3)
    title(main="voom: Mean-variance trend", cex.main=4.5)
    lines(l, col = "red",lwd=3)
    #mtext("Sqrt( standard deviation )", side=2, line=6, cex=4)
    #mtext("log2( count size + 0.5 )", side=1, line=6, cex=4)
    
  }
  f <- approxfun(l, rule = 2)
  if (fit$rank < ncol(design)) {
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[, 
                                                                  j, drop = FALSE])
  }
  else {
    fitted.values <- fit$coef %*% t(fit$design)
  }
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
  fitted.logcount <- log2(fitted.count)
  w <- 1/f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)
  out$E <- y
  out$weights <- w
  out$design <- design
  if (is.null(out$targets)) 
    out$targets <- data.frame(lib.size = lib.size)
  else out$targets$lib.size <- lib.size
  if (save.plot) {
    out$voom.xy <- list(x = sx, y = sy, xlab = "log2( count size + 0.5 )", 
                        ylab = "Sqrt( standard deviation )")
    out$voom.line <- l
  }
  new("EList", out)
}

print("template_function_DEG_VB013.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_VB013_LCF<-readRDS(paste0(rds_output,"/var_VB013_LCF.rds"))
Input_is_Seurat_count <- 0
for(item in var_VB013_LCF){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_VB013_LCF<-as.data.frame(var_VB013_LCF)}else{var_VB013_LCF <- var_VB013_LCF}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Meta_VB013_Outlier_rem<-readRDS(paste0(rds_output,"/var_Meta_VB013_Outlier_rem.rds"))
Input_is_Seurat_count <- 0
for(item in var_Meta_VB013_Outlier_rem){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Meta_VB013_Outlier_rem<-as.data.frame(var_Meta_VB013_Outlier_rem)}else{var_Meta_VB013_Outlier_rem <- var_Meta_VB013_Outlier_rem}
invisible(graphics.off())
var_DEG_VB013<-DEG_VB013(var_VB013_LCF,var_Meta_VB013_Outlier_rem)
invisible(graphics.off())
saveRDS(var_DEG_VB013, paste0(rds_output,"/var_DEG_VB013.rds"))
