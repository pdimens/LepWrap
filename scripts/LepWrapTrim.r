#! /usr/bin/env Rscript

suppressMessages(if (!require("tidyverse")) install.packages("tidyverse"))
suppressMessages(library("tidyverse"))
args <- commandArgs(trailingOnly = TRUE)
#args <- c("~/ordered.1", 10, 15)

# args[1] is the OrderMarkers2 output file
# args[2] is the centimorgan cutoff distance
# args[3] is the % of edge markers to scan

lgfile <- read.delim(
  args[1], 
  header = FALSE, 
  sep = "\t", 
  comment.char="#"
) %>%
  mutate(Mpass = TRUE, Fpass = TRUE)

## setup output file names ##
# split the filename by path
filename <- unlist(strsplit(args[1], "/"))
# pop out just the filename
filename <- filename[length(filename)]
lg <- unlist(strsplit(filename, "\\."))[2]

#========= output instantiation ========#
outfile_base <- paste("5_Trim", filename, sep = "/")
outfile_log_base <- paste("5_Trim", "logs", filename, sep = "/")
plotfile_base <- paste("5_Trim", "plots", filename, sep = "/")
plotfile <- paste(plotfile_base, "trim.pdf", sep = ".")

##### Pruning the ends #####
dist_thresh <- as.numeric(args[2])
# if the percent threshold is given as an interger, convert it to a decimal
edge_length <- as.numeric(args[3])
if(edge_length >= 1){
  edge_length <- edge_length * .01
}

for (j in 2:3){   # iterate over male (2) and female (3)
  # trim beginning
  n_markers <- length(lgfile$V1) * edge_length
  for(a in 1:n_markers){ #first n% of total markers from the beginning
    diff <- abs(lgfile[a+1,j]-lgfile[a,j]) # difference between two points
    if( diff > dist_thresh ){ # is the difference between the two points > distance argument?
      lgfile[1:a, j+4] <- FALSE # mark that marker and all markers BEFORE it as FALSE
    }
  }
  # trim end
  filelen<-length(lgfile$V1)  # get new file lengths for each time we remove NA's
  for(z in filelen:(filelen-n_markers)){  #iterate n% total markers in starting from the end
    diff <- abs(lgfile[z,j]-lgfile[z-1,j]) # difference between two points
    if( diff > dist_thresh ){ # is the difference between the two points > distance argument?
      lgfile[filelen:z,j+4] <- FALSE # mark that marker and all markers AFTER it as FALSE
    }
  }
  
  # create new table of markers passing QC
  cleaned_markers <- lgfile %>% filter(Mpass == TRUE & Fpass == TRUE)
} ###


# isolate bad markers
removed_markers <- (lgfile %>% filter(Mpass == FALSE | Fpass == FALSE))$V1 

# get simple counts
rm_male <- (lgfile %>% filter(Mpass == FALSE & Fpass == TRUE))$V1 %>% length()
rm_female <- (lgfile %>% filter(Fpass == FALSE & Mpass == TRUE))$V1 %>% length()
rm_both <- (lgfile %>% filter(Mpass == FALSE & Fpass == FALSE))$V1 %>% length()

QA <- function(x,y){
  if(x & y) return("pass")
  if(!x & !y) return("both")
  if(!x & y) return("male")
  if(x & !y) return("female")
}

QAfix <- function (x,y){
  if(x == "female" && y == "Female") {
    return("fail")
  } else if (x == "male" && y == "Male") {
    return("fail")
  } else {
    return("pass")
  }
}
  
  
plot_df <- lgfile %>%
              rename(Male = V2, Female = V3) %>%
              mutate(Marker = seq_along(Mpass)) %>%
              rowwise() %>%
              mutate(Fail = QA(Mpass, Fpass)) %>%
              select(Marker,Male, Female, Fail) %>%
              pivot_longer(c(Male, Female), names_to = "Sex", values_to = "Position") %>%
              rowwise() %>%
              mutate(Fail = QAfix(Fail, Sex))

# locations for vertical lines
n_markers <- length(lgfile$V1) 
front_edge <- round(n_markers * edge_length, digits = 0)
back_edge <- round(n_markers - (front_edge), digits = 0)

plot_df %>%
  ggplot(aes(Marker, Position)) +
  geom_point(aes(color = Fail), shape = 19) +
  scale_color_manual(values=c("indianred2", "dodgerblue")) +
  geom_vline(xintercept = front_edge, linetype = "dashed", size = 0.2) +
  geom_vline(xintercept = back_edge, linetype = "dashed", size = 0.2) +
  labs(
    title = paste("Edge Cluster Trimming Results for Linkage Group", lg),
    subtitle = paste0("Markers Failing QC: ", rm_female, " female, ", rm_male, " male, ", rm_both, " both (", rm_female+rm_male+rm_both, " total)" ),
    caption = paste0("edge length: ", edge_length*100, "% of markers, cM threshold: ", dist_thresh),
    x = "Marker Number", 
    y = "Position (cM)", 
    color = "QA Status"
    ) +
  facet_wrap(~Sex)

suppressMessages(ggsave(plotfile, width = 7, height = 4, units = "in"))

# outputting filtered files
num_rm <- length(removed_markers)

writeLines(readLines(args[1], n=3), con = paste(outfile_base, "trimmed", sep = "."))
write.table(
  cleaned_markers[,1:5], 
  file = paste(outfile_base, "trimmed", sep = "."), 
  sep = "\t",
  quote = FALSE, 
  row.names = FALSE,
  col.names = FALSE,
  append=TRUE
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