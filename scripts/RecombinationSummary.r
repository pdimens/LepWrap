#! /usr/bin/env Rscript
# This script will parse all the recombination logs of LepMap3/LepMak3r
# and output a table of summary statistics of recombination of samples
suppressMessages(if (!require("dplyr")) install.packages("dplyr"))
suppressMessages(library("dplyr"))

# setup outfile
args = commandArgs(trailingOnly = TRUE)
# path = ordermarkers or reordermarkers (for reusability)
path = args[1]
# working dir = project_folder/[re]ordermarkers/logs/recombinations
setwd(paste(getwd(),path, "logs/recombinations", sep = "/"))

# generate list of recombination files
files <- list.files(paste(path, "logs", "recombination", sep = "/"))

files <- c("~/ordered.1.1.recombinations", "~/ordered.2.1.recombinations")
recomb_df <- data.frame()

for(i in files){ 
  LG <- unlist(strsplit(i, "\\."))[2]
  recomb <- read.csv(i, skip = 1, header = FALSE, sep="")[,c(2,3,5)]
  recomb$LG <- LG
  recomb_df <- rbind(recomb_df, recomb)
}

recomb_summary <- recomb_df %>% group_by(V2,V3) %>% summarise(min(V5), max(V5), mean(V5), sd(V5))
names(recomb_summary) <- c("family", "sample", "minimum", "maximum", "mean", "standard_deviation")

outfile <- paste(getwd(),"recombination.summary", sep = "/")

write.table(
  recomb_summary,
  file=outfile,
  append=FALSE, 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = TRUE
)