########### Main R file calling other R functions and generating final analysis output####################
wd <- getwd()

## import source files
source ("readandclean.r")
source ("imputeandnormalize.r")
source ("modelfitting.r")
source ("visualization.r")


## Read and clean data
df = readandclean("data/proteinGroups_166_Bioinf.txt")

## Apply filtering for rows vithout valid measurments ( all empty cols)
df.F = filter_valids(df,
                     conditions = c("H0", "Hc"),
                     min_count = c(2, 2),
                     at_least_one = TRUE)


## Remove rows where KEEP is FALSE
df.F = filter(df.F, KEEP)

## Normalize data
df.FN = median_centering(df.F)

## Apply imputation
df.FNI = impute_data(df.FN)


## Limma statistics t-test

# define treatment and control groups for two group comparison, assuming 3 cases and 3 controls
ct <- c("LOG2.H0_1", "LOG2.H0_2", "LOG2.H0_3")
tr <- c("LOG2.Hc_1", "LOG2.Hc_2", "LOG2.Hc_3")


fit = MarayLM(df.FNI,tr,ct)
result = fitlimma(df.FNI,tr,ct)
signproteins = result  %>%  # print significant protein list based on p < 0.005
  filter (p.val < 0.005 & logFC >3)


# print and save file 

print (signproteins)
write.csv(signproteins,paste(wd,"signproteins.csv"),row.names=FALSE)


## Visualize data in Volcano plot (FoldChange~-log10(p.value))

volcanoLIMMA(fit)
volcano(result)  # another viz of sign proteins
MDPlot(fit)
QQPlot(fit)



## enrichments analysis done online via:
#https://reactome.org/PathwayBrowser/#/ 

tinytex::parse_packages()





      