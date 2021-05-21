#! /usr/bin/env Rscript

suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))

args <- commandArgs(trailingOnly = TRUE)
# args[1] = mareydata file
# args[2] = marey midpoint data (optional)

allmaps <- suppressMessages(read_tsv(gzfile(args[1]), col_names = FALSE)) %>%
  select(X3, X2, X5) %>%
  rename(lg = X3, Mb = X2, cM = X5) %>%
  mutate(Mb = Mb/1000000) %>%
  group_by(lg) %>%
  mutate(marker = seq_along(Mb))
  
allmaps %>%
  ggplot(aes(x = Mb, y = cM)) +
  geom_point(size = 0.6, color = "dodgerblue")  +
  labs(
      title = "Marker Positions in Anchored+Oriented Assembly",
      subtitle = "The relative marker positions vs their position in the chromosome/LG",
      x = "Physical Position (Mbp)",
      y = "Genetic Position (cM)"
      ) +
  theme(
        legend.position = "top",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, -5, -5, -5),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(colour = "black", size = 11),
        axis.text.y = element_text(colour = "black", size = 11),
        strip.text = element_text(size = 13),
        axis.title.x = element_text(size = 13, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 13, margin = margin(t = 0, r = 5, b = 0, l = 0))
        ) +
  facet_wrap(~lg, ncol = 4, scales = "free")

outfile <- paste0(dirname(args[1]), "/LepAnchor.midpoint.mareymaps.pdf")
savedims <- length(unique(allmaps$lg)) * 2.5

ggsave(outfile, width = 8.5, height = savedims/4, units = "in")

allmaps %>%
  ggplot(aes(x = marker, y = cM)) +
  geom_point(size = 0.6, alpha = 0.5, color = "dodgerblue")  +
  labs(
      title = "Relative Marker Positions within Linkge Groups",
      subtitle =  "The distance of sequential markers from each other in the linkage maps",
      x = "Marker Number",
      y = "Genetic Position (cM)"
      ) +
  theme(
        legend.position = "top",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, -5, -5, -5),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(colour = "black", size = 11),
        axis.text.y = element_text(colour = "black", size = 11),
        strip.text = element_text(size = 13),
        axis.title.x = element_text(size = 13, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 13, margin = margin(t = 0, r = 5, b = 0, l = 0))
        ) +
  facet_wrap(~lg, ncol = 4, scales = "free_x")

outfile2 <- paste0(dirname(args[1]), "/LepAnchor.midpoint.sequentialmaps.pdf")
ggsave(outfile2,  width = 8.5, height = savedims/5, units = "in", limitsize = FALSE)