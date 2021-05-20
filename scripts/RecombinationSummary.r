#! /usr/bin/env Rscript
# This script will parse all the recombination logs of LepWrap
suppressMessages(if (!require("stringr")) install.packages("stringr"))
suppressMessages(library("stringr"))
## setup outfile
# format trailing arguments for script
args <- commandArgs(trailingOnly = TRUE)

# generate list of recombination files
files <- str_sort(list.files(args[1], pattern = "^order", full.name = TRUE), numeric = TRUE)

## read in recombination logs, add LG info, and append to recomb_df
recomb_df <- read.csv(files[1], skip = 1, header = FALSE, sep="")[,c(2,3,5)]
names(recomb_df) <- c("family", "sample", "LG1")

lg <- 2
for (i in files[2:length(files)]){
  # load in file and pull out recomb info
  .recomb <- read.csv(i, skip = 1, header = FALSE, sep="")[,c(2,3,5)]
  names(.recomb) <- c("family", "sample", paste0("LG",lg))

  # add that as a column to the dataframe
  recomb_df <- merge(recomb_df, .recomb, on = c("family", "sample"), all = TRUE)
  # advance column number iterator
  lg <- lg + 1
}

outfile <- paste(args[1], "recombination.summary", sep = "/")

write.table(
  recomb_df,
  file=outfile,
  append=FALSE, 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = TRUE
)