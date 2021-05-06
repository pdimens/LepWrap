#! /usr/bin/env Rscript

suppressMessages(if (!require("dplyr")) install.packages("dplyr"))
suppressMessages(library("dplyr"))
args = commandArgs(trailingOnly = TRUE)
# args[1] is the path to the working directory
# args[2] is either the ordermarkers or reordermarkers folder
# args[3] is the centimorgan cutoff distance
# args[4] is the % of edge markers to scan

path = args[1]
setwd(path)

lgfile <- read.delim(
            args[2], 
            header = FALSE, 
            sep = "\t", 
            comment.char="#"
          )

# split the filename by path
filename <- unlist(strsplit(args[2], "/"))
# pop out just the filename
filename <- filename[length(filename)]

# chop off the iteration number
filename_trunc <- unlist(strsplit(filename, "\\."))
filename_trunc <- paste(filename_trunc[1], filename_trunc[2], sep = ".")

# instantiate QC columns
lgfile$Mpass <- c(TRUE)
lgfile$Fpass <- c(TRUE)
outfile_base <- paste(path, "4_OrderMarkers", "best.trim", filename_trunc, sep = "/")
outfile_log_base <- paste(path, "4_OrderMarkers", "logs/trimming", filename_trunc, sep = "/")

#========= PDF instantiation ========#
PDFPath <- paste(outfile_log_base, "trim.pdf", sep = ".")
suppressMessages(pdf(file=PDFPath, height = 8.5, width = 11)) 
par(mfrow=(c(1,2))) # create 1x2 plots

##### Pruning the ends #####
dist_thresh <- as.numeric(args[3])
# if the percent threshold is given as an interger, convert it to a decimal
edge_length <- as.numeric(args[4])
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
  
  # diagnostic plots
  par(mar=c(3,4.3,2,1)+0.1)  # reduce the padding somewhat
  if(j == 2){
    plot( 
      x = lgfile[,2], 
      bty="n",
      main = "Male",
      col = "darkseagreen", 
      cex = 1.5,
      ylab = paste(filename, "distance"),
      cex.lab = 1.8,
      xlab = ""
    )
    # only plot if markers were removed
    if(length(which(lgfile$Mpass == FALSE)) != 0){
      points(x = which(lgfile$Mpass == FALSE), y = lgfile$V2[lgfile$Mpass == FALSE], col = "indianred2", pch = 13, cex = 1.8 )   # plot bad markers
    }
  } else {
    plot(
      x = lgfile[,j], 
      bty="n",
      col = "deepskyblue3", 
      cex = 1.5,
      main = "Female",
      ylab = "",
      xlab = ""
    )
    # only plot if markers were removed
    if(length(which(lgfile$Fpass == FALSE)) != 0){
      points(x = which(lgfile$Fpass == FALSE), y = lgfile$V3[lgfile$Fpass == FALSE] , col = "indianred2", pch = 13, cex = 1.8 )
    }
  }
}
graphics.off()
#suppressMessages(dev.off())
# isolate bad markers
removed_markers <- (lgfile %>% filter(Mpass == FALSE | Fpass == FALSE))$V1 

# outputting filtered files

num_rm <- length(removed_markers)

writeLines(readLines(args[2], n=3), con = paste(outfile_base, "trimmed", sep = "."))
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