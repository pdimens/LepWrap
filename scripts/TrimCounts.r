#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
# args[1] = trim.details file
# args[2] = number of LG

if (length(readLines(args[1])) < 1){
    tbl <- 0
    vals <- 0
} else {
    tbl <- read.table(args[1], header = FALSE)
    # split lg filenames to pull out lg number
    vals <- as.numeric(sapply(tbl[,1], function(y){strsplit(y, "\\.")[[1]][2]}))
}
# count the occurence of 1:n_lg
counts <- sapply(1:args[2], function(x){sum(vals == x)})
# create dataframe of counts
out_df <- data.frame(lg = 1:args[2], n_removed = counts)
# print to terminal
print(out_df, row.names = FALSE)