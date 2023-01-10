COX_P181 <- function(MacacaMetadata, Norm_P181) {

library(survival)
library(gdata)

P181 <- Norm_P181
survdata <- MacacaMetadata

outvars <- c("mir",  "study", "timepoint", "RH", "LCI", "UCI", "p", "regr_p", "fdr", "fwer", "fdr.regr", "fwer.regr") # "", 
studies <- c("P181")
timepoints <- list()
# timepoints[[1]] <- c("BASAL", "wk 13")
##timepoints[[1]] <- timepoints[[2]] <- timepoints[[3]] <- c("wk 13") # fixing inner loop, kludge
timepoints[[1]] <- c("wk 13") # P181 and VB013

names(timepoints) <- c("P181") # the name of the study to call the list of timepoints
analysis.text <- "10_15_21"
for(study in studies)
{
  thisstudy <-get(study)
  row.names(thisstudy) <- thisstudy[,1]
  thisstudy <- thisstudy[,-1,]
  thisstudy <- as.data.frame(t(thisstudy))
  names(thisstudy) <- gsub("-", "_", names(thisstudy), fixed = T)
  nmirs <- ncol(thisstudy)
  mirs <- names(thisstudy)
  # browser()
  thisstudy$animal <- substr(row.names(thisstudy), 1, 4)
  coxdata <- merge(survdata, thisstudy, by = "animal")
  coxdata$censor <- ifelse(coxdata$TOA == 12, 0, 1)
  for(timepoint in timepoints[[study]]) #"BASAL",
  {
    # browser()
    # temp
    print(c(thisstudy, timepoint))
    coxoutput <- as.data.frame(matrix(nrow = nmirs, ncol = length(outvars), dimnames = list(NULL, outvars)))
    coxoutput$study <- study
    coxoutput$timepoint <- timepoint
    survmir <- Surv(coxdata$TOA, coxdata$censor)
    for (i in 1:nmirs)
    {
      mir <- mirs[i]
      print(c(i, as.character(mir)))
      coxoutput[i, "mir"] <- mir
      if (dim(table(coxdata[,mir])) < 3) next
      coxdata$thismir <- coxdata[, mir]
      table(coxdata$thismir, coxdata$censor)
      coxrslt <- coxph(survmir ~ thismir, data = coxdata)
      if(is.na(coxrslt$coefficients)) next
      smry <- summary(coxrslt)
      # browser()
      coxoutput[i, "RH"] <- exp(coxrslt$coefficients[1])
      coxoutput[i, "LCI"] <- smry$conf.int["thismir", "lower .95"]
      coxoutput[i, "UCI"] <- smry$conf.int["thismir", "upper .95"]
      coxoutput[i, "p"] <- smry$logtest[3]
      rslt <- lm(coxdata$TOA ~ coxdata$thismir)
      anv <- anova(rslt)
      coxoutput[i, "regr_p"] <- anv["coxdata$thismir", "Pr(>F)"]
      # browser()
    }
    # sigout <- subset(coxoutput, p < 0.05 | regr_p < 0.05)
    coxoutput <- subset(coxoutput, !is.na(p))
    coxoutput$fdr <- p.adjust(coxoutput$p, "fdr")
    coxoutput$fwer <- p.adjust(coxoutput$p, "holm")
    coxoutput$fdr.regr <- p.adjust(coxoutput$regr_p, "fdr")
    coxoutput$fwer.regr <- p.adjust(coxoutput$regr_p, "holm")
    coxoutput <- coxoutput[order(coxoutput$p),]
    sigout <- subset(coxoutput, p < 0.05 | regr_p < 0.05)
 ##   write.csv(sigout, file = paste("sigvars cox and regr 2021", study, timepoint, analysis.text,  ".csv", sep = "_"))
 ##   write.csv(coxoutput, file = paste("allvars cox and regr 2021", study, timepoint, analysis.text,  ".csv", sep = "_"))
  }
}

##write.csv(coxall, file = "VB013 tst 9_18_21.csv")
##wk13vb013 <- subset(coxall, study == "VB013" & timepoint == "wk 13")
##wk13vb013$FDRregr.p <- p.adjust(wk13vb013$regr_p)
##write.csv(wk13vb013[order(wk13vb013$regr_p),], file = "wk13vb013rslts matures 32 at 50%.csv")

return(sigout)
}

print("template_function_COX_P181.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_MacacaMetadata<-readRDS(paste0(rds_output,"/var_MacacaMetadata.rds"))
Input_is_Seurat_count <- 0
for(item in var_MacacaMetadata){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_MacacaMetadata<-as.data.frame(var_MacacaMetadata)}else{var_MacacaMetadata <- var_MacacaMetadata}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Norm_P181<-readRDS(paste0(rds_output,"/var_Norm_P181.rds"))
Input_is_Seurat_count <- 0
for(item in var_Norm_P181){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Norm_P181<-as.data.frame(var_Norm_P181)}else{var_Norm_P181 <- var_Norm_P181}
invisible(graphics.off())
var_COX_P181<-COX_P181(var_MacacaMetadata,var_Norm_P181)
invisible(graphics.off())
saveRDS(var_COX_P181, paste0(rds_output,"/var_COX_P181.rds"))
