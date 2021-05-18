#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("## Print the summary of a map to stdout ##")
  cat("\n[usage]: MapSummary.r <map.file>")
  cat("\n\n## Print summary of all map.* files in a directory to file ##")
  cat("\n[usage]: MapSummary.r <map.directory>")
  cat("\n")
  q()
}

# if args[1] is the target mapfile
if (file_test("-f", args[1])) {
  .dat <- read.table(args[1], header = FALSE, skip = 1)
  # make a table of counts
  .counts <- as.data.frame(table(.dat$V1))
  names(.counts) <- c("lg", "count")
  print(.counts, row.names = FALSE)

# if args[1] is the target directory
} else if (file_test("-d", args[1])) {
library(stringr)
targetdir <- args[1]

# pattern match to find the map files
flist <- list.files(targetdir, full.names = TRUE, pattern = "^map.")
# natural language sort
flist <- str_sort(flist, numeric = TRUE)

for (i in flist){
  # pull the file name from the end
  mapfile <- tail(strsplit(i, "\\/")[[1]], n = 1)
  # pull out the LOD value
  lod<- strsplit(mapfile, "\\.")[[1]][2]
  # read in the file
  .dat <- read.table(i, header = FALSE, skip = 1)
  # make a table of counts
  .counts <- as.data.frame(table(.dat$V1))
  # rename columns
  names(.counts) <- c("LG", paste0("LOD.", lod))
  if (i == flist[1]){
    # if it's the first file, instantiate the counts dataframe
    summtable <- .counts
  } else {
    # otherwise do a matchy-matchy (outer join) to update/bind the dataframe
    summtable <- merge(summtable, .counts, by = "LG", all.x = TRUE, all.y = TRUE)
  }
}

# replace NA's as 0
summtable[is.na(summtable)] <- 0
summtable <- summtable[order(summtable$LG),]

# generate output filenames
out_tmp <-paste0(targetdir, "/all.map.summary")
out_file <- paste0(targetdir, "/all.maps.summary")

write.table(summtable, file = out_tmp, quote = FALSE, row.names = FALSE, col.names = TRUE)
# use GNU column command to make the table fixed-width and remove the tmp file
cmd <- paste("column -t", out_tmp, ">", out_file, "&& rm", out_tmp)
system(cmd)

cat(paste0("Examine the map summary (", out_file, ") and decide on the best map before proceeding"))

} else {
  cat("Error: the argument must be either a mapfile or a directory of mapfiles")
}