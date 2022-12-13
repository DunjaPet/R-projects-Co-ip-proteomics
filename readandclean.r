library(dplyr)

# read and preprocess data file

readandclean <-function (file){
  
  library(dplyr) 

  # Read data file
  raw = read.delim(file, stringsAsFactors = FALSE, colClasses = "character")
  
  
  # filter our _RV and _CON peptoides
  
  to.remove <- grep("CON|REV", raw$Protein.IDs, value=TRUE)
  
  
    df <- raw %>%
    filter(!Protein.IDs %in% to.remove)
  
  
  # keeps only the first instance of Gene 
  df$Gene.names = sub(";.*", "", df$Gene.names)
  
  # remove one hit wonders
  df <- subset(df, df$Peptide.counts..all. != 1)
  
  # Extract names of LFQ.intensity columns
  LFQs = grep("^LFQ.intensity", names(df), value = TRUE)
  
  # Cast as numeric
  df[LFQs] = sapply(df[LFQs], as.numeric)
  
  # Assign column names for log2-transformed data
  LFQ.logs = sub("^LFQ.intensity", "LOG2", LFQs)   # rename intensity columns
  
  # log2 transform of LFG intensity
  
  df[LFQ.logs] = log2(df[LFQs])
  
  # select only columns of interest
  
  df <- select (df,c(Gene.names, LFQ.intensity.H0_1 , LFQ.intensity.H0_2, LFQ.intensity.H0_3,LFQ.intensity.Hc_1, LFQ.intensity.Hc_2, LFQ.intensity.Hc_3, LOG2.H0_1,LOG2.H0_2,LOG2.H0_3,LOG2.Hc_1,LOG2.Hc_2,LOG2.Hc_3 ))
  
  return (df)
}

#red <- readandclean("data/proteinGroups_166_Bioinf.txt")
  
