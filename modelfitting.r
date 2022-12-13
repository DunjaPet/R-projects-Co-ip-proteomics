#############  per http://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html #######################

library(limma)
library(qvalue)

# data preprocessing, load all functions
source("http://www.biostat.jhsph.edu/~kkammers/software/eupa/source.functions.r")

fitmodel<- function (df,tr,ct){

  
  # define design according to syntax of limma package
  design <- model.matrix(~factor(c(2,2,2,1,1,1)))
  colnames(design) <- c("Intercept", "Diff")
  
  res.eb <- eb.fit(df.FNI[, c(tr,ct)], design)
  return (res.eb)
}

# fit LIMMA to get MMaray#####

MarayLM <- function (df,tr,ct){
  design <- model.matrix(~factor(c(2,2,2,1,1,1)))
  fit <- lmFit(df.FNI[, c(tr,ct)],design)
  fit <- eBayes(fit)
  return (fit)
}

#### LIMMA from R documentation #######

fitlimma<- function (df,tr,ct){
  design <- model.matrix(~factor(c(2,2,2,1,1,1)))
  fit <- lmFit(df.FNI[, c(tr,ct)],design)
  fit <- eBayes(fit)

  
  # merge limma results with gene names and log2 from data 
  
  logFC <- fit$coefficients[, 2]
  p.val <- fit$p.value[, 2]
  result <- data.frame(logFC,p.val) %>%
    merge(df.FNI,result, by.x = 0, by.y = 0) %>%
    select (Gene.names,starts_with("LOG2"),logFC,p.val) %>%
    arrange(p.val)
  return(result)
}







