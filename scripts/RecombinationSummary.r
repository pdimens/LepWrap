#! /usr/bin/env Rscript

# This script will parse all the recombination logs of LepWrap
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

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

.recomb_df <- recomb_df %>% select(-family, -sample) 
cnames <- names(.recomb_df)
cnames <- c(c("family", "sample", "mean", "max"), cnames)
.recomb_df <- .recomb_df %>%
  mutate(mean = round(rowMeans(., na.rm = TRUE), digits = 2), max = pmap(., max, na.rm = TRUE))

recomb_df$mean <- .recomb_df$mean
recomb_df$max <- .recomb_df$max
recomb_df <- recomb_df[, cnames]

outfile <- paste(args[1], "recombination.summary", sep = "/")

print(recomb_df, row.names = FALSE, width = 10000)