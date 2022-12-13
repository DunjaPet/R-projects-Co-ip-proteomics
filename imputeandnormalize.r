
#source("readandclean.r")


#df <- readandclean("data/proteinGroups_166_Bioinf.txt")

##remove rows with max 2/3 empty values 

filter_valids = function(df, conditions, min_count, at_least_one = FALSE) {
  # df = data frame containing LOG2 data for filtering and organized by data type
  # conditions = a character vector dictating the grouping
  # min_count = a numeric vector of the same length as "conditions" indicating the minimum 
  #     number of valid values for each condition for retention
  # at_least_one = TRUE means to keep the row if min_count is met for at least one condition
  #     FALSE means min_count must be met across all conditions for retention
  
  log2.names = grep("^LOG2", names(df), value = TRUE)   # Extract LOG2 column names
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, log2.names, value = TRUE, perl = TRUE))
  
  cond.filter = sapply(1:length(cond.names), function(i) {
    df2 = df[cond.names[[i]]]   # Extract columns of interest
    df2 = as.matrix(df2)   # Cast as matrix for the following command
    sums = rowSums(is.finite(df2)) # count the number of valid values for each condition
    sums >= min_count[i]   # Calculates whether min_count requirement is met
  })
  
  if (at_least_one) {
    df$KEEP = apply(cond.filter, 1, any)
  } else {
    df$KEEP = apply(cond.filter, 1, all)
  }
  return(df)  # No rows are omitted, filter rules are listed in the KEEP column
}


## Apply filtering
#df.F = filter_valids(df,
                     #conditions = c("H0", "Hc"),
                     #min_count = c(2, 2),
                     #at_least_one = TRUE)



## Remove rows where KEEP is FALSE
#df.F = filter(df.F, KEEP)

## Data normalization function
median_centering = function(df) {
  # df = data frame containing LOG2 columns for normalization
  LOG2.names = grep("^LOG2", names(df), value = TRUE)
  
  df[, LOG2.names] = lapply(LOG2.names, 
                            function(x) {
                              LOG2 = df[[x]]
                              LOG2[!is.finite(LOG2)] = NA   # Exclude missing values from median calculation
                              gMedian = median(LOG2, na.rm = TRUE)
                              LOG2 - gMedian
                            }
  )
  
  return(df)
}


## Normalize data
#df.FN = median_centering(df.F)

## Data imputation function
impute_data = function(df, width = 0.3, downshift = 1.8) {
  # df = data frame containing filtered 
  # Assumes missing data (in df) follows a narrowed and downshifted normal distribution
  
  LOG2.names = grep("^LOG2", names(df), value = TRUE)  # find names starting with "LOG2"
  impute.names = sub("^LOG2", "impute", LOG2.names)    # create new names from LOG2 names with impute
  
  # Create new columns indicating whether the values are imputed
  df[impute.names] = lapply(LOG2.names, function(x) !is.finite(df[, x]))
  
  # Imputation
  set.seed(1)
  df[LOG2.names] = lapply(LOG2.names,
                          function(x) {
                            temp = df[[x]]
                            temp[!is.finite(temp)] = NA
                            
                            temp.sd = width * sd(temp[df$KEEP], na.rm = TRUE)   # shrink sd width
                            temp.mean = mean(temp[df$KEEP], na.rm = TRUE) - 
                              downshift * sd(temp[df$KEEP], na.rm = TRUE)   # shift mean of imputed values
                            
                            n.missing = sum(is.na(temp))
                            temp[is.na(temp)] = rnorm(n.missing, mean = temp.mean, sd = temp.sd)                          
                            return(temp)
                          })
  return(df)
}

## Apply imputation
#df.FNI = impute_data(df.FN)
