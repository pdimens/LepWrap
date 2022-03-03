#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
# args[1] is the OrderMarkers2/LepAnchor output file
# args[2] is the centiMorgan cutoff threshold
# args[3] is the % of edge markers to scan
# args[4] is the output directory

#==== functions to translate QC flags ====#
QA <- function(x,y){
  if(x & y) return("pass")
  if(!x & !y) return("both")
  if(!x & y) return("male")
  if(x & !y) return("female")
}

QAfix <- function (x,y){
  if(x == "female" & y == "Female") {
    return(F)
  } else if (x == "male" & y == "Male") {
    return(F)
  } else if(x == "both") {
    return(F)
  } else {
    return(T)
  }
}

#==== read in file =====#

lgfile <- read.delim(
  args[1],
  header = FALSE, 
  sep = "\t", 
  comment.char="#"
) %>%
  mutate(Mpass = T, Fpass = T)

## setup output file names ##
# split the filename by path
filename <- unlist(strsplit(args[1], "/"))
# pop out just the filename
filename <- filename[length(filename)]

# Modify which columns to look at based on input format
# LM3 input has 5 columns + 2 from M|Fpass
# LA input has 6 columns + 2 from M|Fpass
if(ncol(lgfile) == 7) {
  # LepMap3 input file
  idxcol <- c(2,3)
  keepcol <- 1:5
  lg <- unlist(strsplit(filename, "\\."))[2]
  colshift <- 4
} else {
  # LepAnchor input file
  idxcol <- c(5, 6)
  keepcol <- 1:6
  lg <- unlist(strsplit(filename, "\\."))[3]
  colshift <- 2
}

#========= instantiate output ========#
dir.create(args[4], showWarnings = FALSE)
dir.create(paste0(args[4],"/plots"), showWarnings = FALSE)
dir.create(paste0(args[4],"/logs"), showWarnings = FALSE)
dir.create(paste0(args[4],"/QC_raw"), showWarnings = FALSE)
outfile_base <- paste(args[4], filename, sep = "/")
outfile_log_base <- paste(args[4], "logs", filename, sep = "/")
plotfile_base <- paste(args[4], "plots", filename, sep = "/")
plotfile <- paste(plotfile_base, "trim.pdf", sep = ".")
rawfile_base <- paste(args[4], "QC_raw", filename, sep = "/")


#====== Pruning the ends ======#
# if the percent threshold is given as an integer, convert it to a decimal
dist_thresh <- as.numeric(args[2])
if(dist_thresh >= 1){
  dist_thresh <- dist_thresh * .01
}

dist_thresh_all <- c(
  abs(max(lgfile[, idxcol[1]]) - min(lgfile[, idxcol[1]])) * dist_thresh,     # male
  abs(max(lgfile[, idxcol[2]]) - min(lgfile[, idxcol[2]])) * dist_thresh      # female
)

edge_length <- as.numeric(args[3])
if(edge_length >= 1){
  edge_length <- edge_length * .01
}

n_markers <- nrow(lgfile)
forward_start <- round(n_markers * edge_length, digits = 0)
reverse_start <- round(n_markers - forward_start, digits = 0)

# iterate over male (2) and female (3)
threshidx <- 1
for (j in idxcol){
  # sort on column
  lgfile <- arrange(lgfile, j)
  dist_thresh <- dist_thresh_all[threshidx]
  # trim beginning
  for(a in forward_start:2){ #first n% of total markers from the beginning
    diff <- abs(lgfile[a,j]-lgfile[a-1,j]) # difference between two points
    if( diff > dist_thresh ){ # is the difference between the two points > distance argument?
      lgfile[(a-1):1, j + colshift] <- FALSE # all markers BEFORE it as FAIL
      break()
    }
  }
  # trim end
  for(z in reverse_start:(n_markers-1)){ #last n% total markers starting from the back edge going out
    diff <- abs(lgfile[z+1,j]-lgfile[z,j]) # difference between two points
    if( diff > dist_thresh ){ # is the difference between the two points > distance argument?
      lgfile[(z+1):n_markers, j + colshift] <- FALSE # all markers AFTER it as FAIL
      break()
    }
  }
  threshidx <- threshidx + 1
}

# isolate bad markers
removed_markers <- lgfile %>% filter(!Mpass | !Fpass) %>% select(1) 

# get simple counts
rm_male <- lgfile %>% filter(!Mpass & Fpass) %>%  nrow()
rm_female <- lgfile %>% filter(!Fpass & Mpass) %>% nrow()
rm_both <- lgfile %>% filter(!Mpass & !Fpass) %>% nrow()


#====== Plotting ======#
pdf(NULL)

plot_df <- lgfile %>%
  rename(Male = idxcol[1], Female = idxcol[2]) %>%
  arrange(Male) %>%
  mutate(Marker = seq_along(Mpass)) %>%
  rowwise() %>%
  mutate(Fail = QA(Mpass, Fpass)) %>% 
  select(Marker,Male, Female, Fail) %>%
  pivot_longer(c(Male, Female), names_to = "Sex", values_to = "Position") %>%
  rowwise() %>%
  mutate(Fail = QAfix(Fail, Sex))

plot_df %>%
  ggplot(aes(Marker, Position, color = Fail, shape = Fail)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("dodgerblue", "indianred2"), limits = c(T, F)) +
  scale_shape_manual(values = c(19, 17), limits = c(T, F), guide = "none") +
  geom_vline(xintercept = forward_start, linetype = "dashed", size = 0.2) +
  geom_vline(xintercept = reverse_start, linetype = "dashed", size = 0.2) +
  guides(color = guide_legend(override.aes = list(size = c(2,2), shape = c(19, 17)))) +
  labs(
    title = paste("Edge Cluster Trimming Results for Linkage Group", lg),
    subtitle = paste0("Markers Failing QC: ", rm_female, " female, ", rm_male, " male, ", rm_both, " both (", rm_female+rm_male+rm_both, " total)" ),
    caption = paste0(edge_length*100, "% of edge markers, ", args[2], "% cM threshold: ", round(dist_thresh_all[2], digits = 2), "(F)  &  ", round(dist_thresh_all[1], digits = 2), "(M)"),
    x = "Marker Number", 
    y = "Position (cM)", 
    color = "Pass QA"
    ) +
  facet_wrap(~Sex)

suppressMessages(ggsave(plotfile, width = 7, height = 4, units = "in"))

#====== outputting filtered files ======#
# get an overall pass/fail for each marker
lgfile$QC <- lgfile$Mpass & lgfile$Fpass
# get the indices of the fails
fail_idx <- which(!lgfile$QC)
# prepend a comment to the flagged markers so LepMap3 ignores them
lgfile[fail_idx, 1] <- paste0("#", lgfile[fail_idx, 1])
# re-scale remaining markers to 0 by subtracting the minimum genetic position for each sex
lgfile[,idxcol[1]] <- round(lgfile[,idxcol[1]] - (min(lgfile[lgfile$QC,idxcol[1]])), digits = 3)
lgfile[,idxcol[2]] <- round(lgfile[,idxcol[2]] - (min(lgfile[lgfile$QC,idxcol[2]])), digits = 3)

# write header to a new file
writeLines(readLines(args[1], n=2), con = paste(outfile_base, "trimmed", sep = "."))
# write the remainder to that file
write.table(
  lgfile[, keepcol], 
  file =  paste(outfile_base, "trimmed", sep = "."), 
  sep = "\t",
  quote = FALSE, 
  row.names = FALSE,
  col.names = FALSE,
  append=TRUE
)

write.table(
  lgfile, 
  file = paste(rawfile_base, "filtered.raw", sep = "."), 
  sep = "\t",
  quote = FALSE, 
  row.names = FALSE,
  col.names = FALSE,
)

write.table(
  removed_markers,
  file=paste(outfile_log_base, "removed", sep = "."),
  append=FALSE, 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = FALSE
)