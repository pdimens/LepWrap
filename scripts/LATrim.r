#! /usr/bin/env Rscript

suppressMessages(if (!require("tidyverse")) install.packages("tidyverse"))
suppressMessages(library("tidyverse"))
suppressMessages(if (!require("cowplot")) install.packages("cowplot"))
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)
# args[1] is the OrderMarkers2 output file
# args[2] is the centimorgan cutoff distance
# args[3] is the % of edge markers to scan
# args[4] is the name of the output folder

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
lg <- (strsplit(lgfile$V1[1], "LG") %>% unlist())[2] %>% as.numeric()

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

##### Pruning the ends #####
dist_thresh <- as.numeric(args[2])
if(dist_thresh >= 1){
  dist_thresh <- dist_thresh * .01
}
.dist_thresh <- dist_thresh * 100
dist_thresh_m <- round(abs(max(lgfile$V5) - min(lgfile$V5)) * dist_thresh, digits = 2)
dist_thresh_f <- round(abs(max(lgfile$V6) - min(lgfile$V6)) * dist_thresh, digits = 2)


# if the percent threshold is given as an integer, convert it to a decimal
edge_length <- as.numeric(args[3])
if(edge_length >= 1){
  edge_length <- edge_length * .01
}
.edge_length <- edge_length * 100
n_markers <- length(lgfile$V1)

forward_start <- round(n_markers * edge_length, digits = 0)
reverse_start <- round(n_markers - forward_start, digits = 0)

for (j in 5:6){   # iterate over male (5) and female (6)
  # sort on column
  if (j == 5){
    lgfile <- arrange(lgfile, V5)
    dist_thresh <- dist_thresh_m
  } else {
    lgfile <- arrange(lgfile, V6)
    dist_thresh <- dist_thresh_f
  }
  # trim beginning
  for(a in forward_start:2){ #first n% of total markers starting from the front edge, going out
    diff <- abs(lgfile[a,j]-lgfile[a-1,j]) # difference between two points
    if( diff > dist_thresh ){ # is the difference between the two points > distance argument?
      lgfile[(a-1):1, j+2] <- F # mark that marker and all markers BEFORE it as FAIL
      break()
    }
  }
  # trim end
  for(z in reverse_start:(n_markers-1)){ #last n% total markers starting from the back edge going out
    diff <- abs(lgfile[z+1,j]-lgfile[z,j]) # difference between two points
    if( diff > dist_thresh ){ # is the difference between the two points > distance argument?
      lgfile[(z+1):n_markers,j+2] <- F # mark that marker and all markers AFTER it as FAIL
      break()
    }
  }
}

# create new table of markers passing QC
cleaned_markers <- (lgfile %>% filter(Mpass & Fpass))[,1:6]
# re-scale cleaned markers to 0 by subtracting the minimum genetic position
cleaned_markers <- cleaned_markers %>%
                    mutate(V5 = V5 - min(V5), V6 = V6 - min(V6))

# isolate bad markers
removed_markers <- (lgfile %>% filter(!Mpass | !Fpass))[,1:6]

# get simple counts
rm_male <- lgfile %>% filter(!Mpass & Fpass) %>% nrow()
rm_female <- lgfile %>% filter(!Fpass & Mpass) %>% nrow()
rm_both <- lgfile %>% filter(!Mpass & !Fpass) %>% nrow()

pdf(NULL)

plot_male <- lgfile %>% arrange(V5) %>%
  ggplot(aes(x = seq_along(V5), y = V5, color = Mpass)) +
  geom_point(shape = 19) +
  scale_color_manual(values = c("dodgerblue", "indianred2"), limits = c(T, F)) +
  geom_vline(xintercept = forward_start, linetype = "dashed", size = 0.2) +
  geom_vline(xintercept = reverse_start, linetype = "dashed", size = 0.2) +
  labs(
    title = "",
    subtitle = paste0(rm_male, " male markers >", dist_thresh_m, "cM trimmed"),
    caption = paste0(.edge_length, "% edge markers, ", .dist_thresh, "% cM"),
    x = "Marker Number", 
    y = "Position (cM)", 
    color = "Pass Filtering"
  )

plot_female <- lgfile %>% arrange(V6) %>%
  ggplot(aes(x = seq_along(V6), y = V6, color = Fpass)) +
  geom_point(shape = 19) +
  scale_color_manual(values = c("dodgerblue", "indianred2"), limits = c(T, F)) +
  geom_vline(xintercept = forward_start, linetype = "dashed", size = 0.2) +
  geom_vline(xintercept = reverse_start, linetype = "dashed", size = 0.2) +
  labs(
    title = paste("Edge Cluster Trimming for LG:", lg),
    subtitle = paste0(rm_female, " female markers >", dist_thresh_f, "cM trimmed"),
    caption = paste0("Markers failing both M+F: ", rm_both),
    x = "Marker Number", 
    y = "Position (cM)", 
    color = "Pass Filtering",
    legend.position = "none"
  )

plot_grid(plot_female, plot_male, ncol = 2, nrow = 1)

suppressMessages(ggsave(plotfile, width = 7, height = 3, units = "in"))

write.table(
  cleaned_markers, 
  file = paste(outfile_base, "trimmed", sep = "."), 
  sep = "\t",
  quote = FALSE, 
  row.names = FALSE,
  col.names = FALSE,
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