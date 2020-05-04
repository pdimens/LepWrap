#! /usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

suppressMessages(if (!require("dplyr")) install.packages("dplyr"))
suppressMessages(library("dplyr"))
path = args[1]
setwd(args[1])
file.names <- scan(paste(path, args[2], sep = "/"), character(), quote = "")
#file.names <- file.names[order(nchar(file.names), file.names)] #sort by LG
PDFPath <- paste(path, "/best.trimmed/trimming.plots.pdf", sep = "/")
pdf(file=PDFPath, height = 11, width = 8.5) 
par(mfrow=(c(4,2))) # create 4x2 plots

##### Pruning the ends #####
for(i in file.names){  
  lgfile <- read.delim(
    paste(path, i, sep = "/"), 
    header = FALSE, 
    sep = "\t", 
    comment.char="#"
  )
  
  # instantiate QC columns
  lgfile$Mpass <- c(TRUE)
  lgfile$Fpass <- c(TRUE)
  
  for (j in 2:3){   # iterate over male (2) and female (3)
    # trim beginning
    filelength15 <- length(lgfile$V1) * 0.15
    for(a in 1:filelength15){ #first 15% of total markers from the beginning
      diff <- abs(lgfile[a+1,j]-lgfile[a,j]) # difference between two points
      if( diff > 10 ){ # is the difference between the two points > distance argument?
        lgfile[1:a, j+4] <- FALSE
      }
    }
    # trim end
    filelen<-length(lgfile$V1)  # get new file lengths for each time we remove NA's
    for(z in filelen:(filelen-filelength15)){  #iterate 15% total markers in starting from the end
      diff <- abs(lgfile[z,j]-lgfile[z-1,j]) # difference between two points
      if( diff > 10 ){ # is the difference between the two points > distance argument?
        lgfile[filelen:z,j+4] <- FALSE # mark that marker and all markers AFTER it as NA
      }
    }
    
    # create new table of markers passing QC
    cleaned_markers <- lgfile %>% filter(Mpass == TRUE & Fpass == TRUE)  
    
    # diagnostic plots
    par(mar=c(3,4.3,2,1)+0.1)  # reduce the padding somewhat
    just_ordernum <- tools::file_path_sans_ext(i)  # remove the extension from files for plots
    if( j==2 ){
      plot( x = lgfile[,2], 
            bty="n",
            main = "Male",
            col = "slategray", 
            cex = 1.5,
            ylab = paste(just_ordernum, "distance"),
            cex.lab = 1.8,
            xlab = "",
      )
      points(x = which(lgfile$Mpass == FALSE), y = lgfile$V2[lgfile$Mpass == FALSE], col = "coral3", pch = 19, cex = .8 )   # plot bad markers
    } else {
      plot( x = lgfile[,j], 
            bty="n",
            col = "slategray", 
            cex = 1.5,
            main = "Female",
            ylab = "",
            xlab = ""
      )
      points(x = which(lgfile$Fpass == FALSE), y = lgfile$V3[lgfile$Fpass == FALSE] , col = "coral3", pch = 19, cex = .8 )
    }
  }
  
  # isolate bad markers
  removed_markers <- (lgfile %>% filter(Mpass == FALSE | Fpass == FALSE))$V1 
  
  # outputting filtered files
  filename<- paste("trimmed",i, sep=".")
  print(paste("Removing",length(removed_markers),"markers from",i , "and writing new file", filename, sep = " "))
  writeLines(readLines(i, n=3),con = filename)
  write.table(cleaned_markers[,1:5], 
              file = filename, 
              sep = "\t",
              quote = FALSE, 
              row.names = FALSE,
              col.names = FALSE,
              append=TRUE
  )
  write.table(removed_markers,
              file="bad.markers.txt",
              append=TRUE, 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE, 
              col.names = FALSE
  )
}
suppressMessages(dev.off())
