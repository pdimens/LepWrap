#! /usr/bin/env Rscript

library(stringr)

args <- commandArgs(trailingOnly = TRUE)

# args[1] is the target directory

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

# generate output filenames
out_tmp <-paste0(targetdir, "/all.map.summary")
out_file <- paste0(targetdir, "/all.maps.summary")

write.table(summtable, file = out_tmp, quote = FALSE, row.names = FALSE, col.names = TRUE)
# use GNU column command to make the table fixed-width and remove the tmp file
cmd <- paste("column -t", out_tmp, ">", out_file, "&& rm", out_tmp)
system(cmd)

cat(paste0("Examine the map summary (", out_file, ") and decide on the best map before proceeding"))
cat("\nIf using a screen/tmux environment, detach this session and return to it with the appropriate command when ready\n")

### these plots are completely unnecessary, but I'll leave them here anyway
#plot_table <- summtable %>% 
#              pivot_longer(cols = !LG, names_to = "LOD", values_to = "n_markers") %>% 
#              filter(n_markers > 0) %>%
#              mutate(LOD = as.numeric(LOD))

#plot_table %>% ggplot(aes(LOD, LG,  label = n_markers)) +
#               geom_text() +
#               scale_x_continuous("LOD", labels = as.character(plot_table$LOD), breaks = plot_table$LOD)
