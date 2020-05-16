#! /usr/bin/env Rscript
# This script will parse all the recombination logs of LepMap3/LepMak3r
# and output a table of summary statistics of recombination of samples
suppressMessages(if (!require("dplyr")) install.packages("dplyr"))
suppressMessages(library("dplyr"))

## setup outfile
# format trailing arguments for script
args = commandArgs(trailingOnly = TRUE)
# path = ordermarkers or reordermarkers (for reusability)
path = args[1]
original_wd <- getwd()
# working dir = project_folder/[re]ordermarkers/logs/recombinations
setwd(paste(getwd(),path, "logs/recombination", sep = "/"))

# generate list of recombination files
files <- list.files(getwd())

# instantiate empty dataframe that will hold all the recombination info
recomb_df <- data.frame()

## read in recombination logs, add LG info, and append to recomb_df
for(i in files){
  # pull out linkage group number from filename
  LG <- unlist(strsplit(i, "\\."))[2]
  # read in file and keep only family|sample|recombinations
  recomb <- read.csv(i, skip = 1, header = FALSE, sep="")[,c(2,3,5)]
  # add column identifying the linkage group
  recomb$LG <- LG
  # append the formatted dataframe to the total dataframe
  recomb_df <- rbind(recomb_df, recomb)
}

## summary information
# calculate min, max, mean, and sd across iterations for each individual in each LG
recomb_summary <- recomb_df %>% group_by(V2,V3) %>% summarise(min(V5), max(V5), format(round(mean(V5), digits = 4), nsmall = 4), format(round(sd(V5), digits = 4), nsmall = 4)

# rename columns for clarity
names(recomb_summary) <- c("family", "sample", "minimum", "maximum", "mean", "standard_deviation")

# generate output file's name based on whether it's re/ordermarkers
outfile <- paste(original_wd, path, "recombination.summary", sep = "/")

## write summary to file
write.table(
  recomb_summary,
  file=outfile,
  append=FALSE, 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = TRUE
)