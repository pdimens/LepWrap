#! /usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
# args[1] is the summary file

# read in the summary table and pop out the LG number
tbl <- read_table(args[1]) %>% suppressMessages()
pdf(NULL)

tbl %>%
  ggplot(aes(lg, n_removed)) +
  geom_segment(
    aes(lg, n_removed, xend = lg, yend=0),
    size = 5,
    color = "#78bcff"
    ) +
  geom_point(shape = 21, fill = "#4aa6ff", color = "transparent", size = 5) +
  scale_x_continuous(labels = as.character(tbl$lg), breaks = tbl$lg) +
  labs(
    title = "Number of markers removed via trimming",
    x = "Linkage Group",
    y = "Number removed"
       ) +
    scale_y_continuous(limits =   c(-0.05, max(tbl$n_removed)+10))

ggsave(paste0(args[1], ".pdf"), height = 4, width = 7, units = "in")
ggsave(paste0(args[1], ".svg"), height = 4, width = 7, units = "in")