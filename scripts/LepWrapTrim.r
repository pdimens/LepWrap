#! /usr/bin/env Rscript

suppressMessages(if (!require("tidyverse")) install.packages("tidyverse"))
suppressMessages(library("tidyverse"))

args <- commandArgs(trailingOnly = TRUE)
# args[1] is the OrderMarkers2 output file
# args[2] is the centiMorgan cutoff threshold
# args[3] is the % of edge markers to scan
# args[4] is the output directory

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
lg <- unlist(strsplit(filename, "\\."))[2]

#========= output instantiation ========#
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
  abs(max(lgfile[, 2]) - min(lgfile[, 2])) * dist_thresh,     # male
  abs(max(lgfile[, 3]) - min(lgfile[, 3])) * dist_thresh      # female
)

edge_length <- as.numeric(args[3])
if(edge_length >= 1){
  edge_length <- edge_length * .01
}

n_markers <- length(lgfile$V1)
forward_start <- round(n_markers * edge_length, digits = 0)
reverse_start <- round(n_markers - forward_start, digits = 0)

for (j in 2:3){   # iterate over male (2) and female (3)
  # sort on column
  lgfile <- arrange(lgfile, j)
  dist_thresh <- dist_thresh_all[j-1]
  # trim beginning
  for(a in forward_start:2){ #first n% of total markers from the beginning
    diff <- abs(lgfile[a,j]-lgfile[a-1,j]) # difference between two points
    if( diff > dist_thresh ){ # is the difference between the two points > distance argument?
      lgfile[(a-1):1, j+4] <- FALSE # all markers BEFORE it as FAIL
      break()
    }
  }
  # trim end
  for(z in reverse_start:(n_markers-1)){ #last n% total markers starting from the back edge going out
    diff <- abs(lgfile[z+1,j]-lgfile[z,j]) # difference between two points
    if( diff > dist_thresh ){ # is the difference between the two points > distance argument?
      lgfile[(z+1):n_markers,j+4] <- FALSE # all markers AFTER it as FAIL
      break()
    }
  }
}

# isolate bad markers
removed_markers <- lgfile %>% filter(!Mpass | !Fpass) %>% select(1) 

# get simple counts
rm_male <- lgfile %>% filter(!Mpass & Fpass) %>%  nrow()
rm_female <- lgfile %>% filter(!Fpass & Mpass) %>% nrow()
rm_both <- lgfile %>% filter(!Mpass & !Fpass) %>% nrow()

# functions to translate QC flags for plotting and other information
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
  } else if(x == "both" & (y == "Male" | y == "Female")) {
    return(F)
  } else {
    return(T)
  }
}

#====== Plotting ======#
pdf(NULL)

plot_df <- lgfile %>%
  rename(Male = V2, Female = V3) %>%
  arrange(Male) %>%
  mutate(Marker = seq_along(Mpass)) %>%
  rowwise() %>%
  mutate(Fail = QA(Mpass, Fpass)) %>% 
  select(Marker,Male, Female, Fail) %>%
  pivot_longer(c(Male, Female), names_to = "Sex", values_to = "Position") %>%
  rowwise() %>%
  mutate(Fail = QAfix(Fail, Sex))

plot_df %>%
  ggplot(aes(Marker, Position, color = Fail)) +
  geom_point(shape = 19, alpha = 0.3) +
  scale_color_manual(values = c("dodgerblue", "indianred2"), limits = c(T, F)) +
  geom_vline(xintercept = forward_start, linetype = "dashed", size = 0.2) +
  geom_vline(xintercept = reverse_start, linetype = "dashed", size = 0.2) +
  labs(
    title = paste("Edge Cluster Trimming Results for Linkage Group", lg),
    subtitle = paste0("Markers Failing QC: ", rm_female, " female, ", rm_male, " male, ", rm_both, " both (", rm_female+rm_male+rm_both, " total)" ),
    caption = paste0(edge_length*100, "% of edge markers, ", args[2], "% cM threshold: ", dist_thresh_all[1], "(M) & ", dist_thresh_all[2], "(F)"),
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
flag_idx <- which(!lgfile$QC)
# prepend a comment to the flagged markers so LepMap3 ignores them
lgfile[flag_idx, 1] <- paste0("#", lgfile[flag_idx, 1])

# write header to a new file
writeLines(readLines(args[1], n=3), con = paste(outfile_base, "trimmed", sep = "."))
# write the remainder to that file
write.table(
  lgfile[,1:5], 
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