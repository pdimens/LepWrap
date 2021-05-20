#! /usr/bin/env Rscript

suppressMessages(library(tidyverse, quietly = TRUE))

args <- commandArgs(trailingOnly = TRUE)
# args[1] is the summary file
# args[2] is the expected number of linkage groups (to fill in zeroes)

# if the summary file is empty (b/c nothing was trimmed)
if (length(readLines(args[1])) < 2){
    tbl <- data.frame(map = 1:args[2], n_removed = 0)
} else {
    # read in the summary table and pop out the LG number
    tbl <- read_table(args[1]) %>%
            rowwise() %>%
            mutate(map = as.numeric(strsplit(map, "\\.")[[1]][2])) %>%
            suppressMessages()

    # find the max number of linkage groups
    lg_num <- 1:args[2]
    missing_lg <- lg_num[!(lg_num %in% tbl$map)]
    missings <- data.frame(n_removed = 0, map = missing_lg)
    tbl <- merge(tbl, missings, on = map, all = TRUE)
}

tbl %>%
  ggplot(aes(map, n_removed)) +
  geom_segment(
    aes(map, n_removed, xend = map, yend=0),
    size = 5,
    color = "#78bcff"
    ) +
  geom_point(shape = 21, fill = "#4aa6ff", color = "transparent", size = 5) +
  scale_x_continuous(labels = as.character(tbl$map), breaks = tbl$map) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title = "Number of markers removed via trimming",
    x = "Linkage Group",
    y = "Number removed"
       )
  
ggsave(paste0(args[1], ".pdf"), height = 4, width = 7, units = "in")
ggsave(paste0(args[1], ".svg"), height = 4, width = 7, units = "in")