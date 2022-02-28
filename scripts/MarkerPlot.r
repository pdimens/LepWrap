#! /usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

pdf(NULL)

lgfile <- read.delim(
  args[1], 
  header = FALSE, 
  sep = "\t", 
  comment.char="#"
) %>%
  select(V2,V3) %>%
  rename(male = V2, female = V3) %>%
  arrange(male) %>% 
  pivot_longer(c(male, female), names_to = "sex", values_to = "position") %>%
  group_by(sex) %>%
  mutate(marker = seq_along(position))

lgfile %>%
  ggplot(aes(marker, position)) +
  geom_point(shape=20, alpha = 0.5, aes(color = sex)) +
  theme(legend.position = "none") +
  facet_wrap(~sex)

outfile <- paste0(args[1], ".seqplot.png")
cat(paste0("Saving plot to ", outfile))
cat("\n")

suppressMessages(ggsave(outfile))