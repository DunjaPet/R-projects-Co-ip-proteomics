
## Viz from LIMMA pckg

volcanoLIMMA <- function(fit){
  volcanoplot(fit,coef=2,highlight=2)
}

# Mean-difference plot
MDPlot<- function (fit){
  plotMD(fit,column=2)
}

# Q-Q plot of moderated t-statistics
QQPlot <- function (fit) {
  qqt(fit$t[,2],df=fit$df.residual+fit$df.prior)
abline(0,1)
}



# volcano ggplot

library(ggplot2)

volcano <- function (result){
p <- ggplot(data = result,
            mapping = aes(x = logFC, y = -log10(p.val), colour= factor(p.val<0.005 & logFC >3) ))

p + geom_point()+
  geom_text(data = result[1:2,],
            mapping = aes(label = Gene.names),hjust=0, vjust=0.5, colour="red",size=3) +
  geom_hline(yintercept = -log10(0.005), size = 1.4, color = "gray80") +
  geom_vline(xintercept = 3, size = 1.4, color = "gray80") +
  geom_vline(xintercept = -3, size = 1.4, color = "gray80")
  
  #ggtitle("Volcano Plot")
  
}



