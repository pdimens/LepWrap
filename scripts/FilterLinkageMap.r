# This R file performs an adaptive method of filtering a linkage map
# It works by creating a spline on the linkage map, then performing
# a sliding window analysis on the residuals of the spline, calculating
# a local threshold of 2 x mean(local residuals) to flag markers as
# possible outliers. 

# More specifically, lgFilter creates two splines ("higher" and 
# "lower" resolution), and performs a 50% overlapping sliding
# window analysis with a different window size for each resolution.
# Pass/Fails go into a simple penalty matrix and markers >= 3 penalties are
# finally flagged as FAIL. 3 fails means the marker failed both times for one
# sliding window and at least once for the other. End-markers that only get
# tested 2x (lack of overlaps) incur double-penalties.

library(tidyverse)
library(splines)

# will require un-gzipped mareydata
marey <- read_delim("data.marey.trimmed", delim = "\t", col_names = FALSE) %>%
  select(X3, X2, X5, X6) %>%
  rename(LG = X3, MB = X2, MALE = X5, FEMALE = X6) %>%
  rowwise() %>%
  mutate(AVG = mean(c(MALE,FEMALE))) %>%
  pivot_longer(c(MALE,FEMALE), values_to = "CM", names_to = "SEX")

lgFilter <- function(mb,cm){
  n <- length(mb)
  knots <- round(n/100)
  .spline <- lm(cm ~ bs(mb, knots = knots))
  res <- .spline$residuals
  # mb = physical position vector 
  # cm = genetic position vector
  # divide by 10 and round to next even integer
  windowsize <- 2 * round((n/10)/2)
  slide <- windowsize/2
  # sliding window of `windowsize`, sliding at half its length each iteration
  # the input vector is divided into 10 ~even windows
  windowseq <- seq(from = 1, to = n-slide, by = slide )
  # create empty penalty matrix
  penalty <- rep(0, n)
  for (.window in windowseq){
    # workaround for last window not being an even size
    if (.window == windowseq[length(windowseq)]){
      .residuals <- res[.window:n]
    }  else{ 
      .residuals <- res[.window:(.window + windowsize-1)]
    }
    threshold <- mean(abs(.residuals)) * 2
    # workaround like above
    # add +1 penalty to end-markers b/c they dont get scanned 4x like the others
    if (.window == windowseq[length(windowseq)]){
      penalty[.window:n] <- penalty[.window:n] + (abs(.residuals) > threshold)
      #slide window up to the halfway point
      halfwindow = median(.window:n)
      # double the penalty for the end points that only get scanned 2x rather than 4x
      penalty[halfwindow:n] <- penalty[halfwindow:n] + 1
    } else {
      penalty[.window:(.window+windowsize-1)] <- penalty[.window:(.window+windowsize-1)] + (abs(.residuals) > threshold)
    }
  }
  # a second window ~30% smaller
  knots <- round(n/70)
  .spline <- lm(cm ~ bs(mb, knots = knots))
  res <- .spline$residuals
  windowsize <- 2 * round((n/7)/2)
  slide <- windowsize/2
  windowseq <- seq(from = 1, to = n-slide, by = slide )
  for (.window in windowseq){
    if (.window == windowseq[length(windowseq)]){
      .residuals <- res[.window:n]
    }  else{ 
      .residuals <- res[.window:(.window + windowsize-1)]
    }
    threshold <- mean(abs(.residuals)) * 2
    # workaround like above
    if (.window == windowseq[length(windowseq)]){
      penalty[.window:n] <- penalty[.window:n] + (abs(.residuals) > threshold)
      #slide window up to the halfway point
      halfwindow = median(.window:n)
      # double the penalty for the end points that only get scanned 2x rather than 4x
      penalty[halfwindow:n] <- penalty[halfwindow:n] + 1
    } else {
      penalty[.window:(.window+windowsize-1)] <- penalty[.window:(.window+windowsize-1)] + (abs(.residuals) > threshold)
    }
  }
  qual <- penalty >= 3
  # convert to PASS/FAIL instead of TRUE/FALSE
  qual[qual] <- "FAIL"
  qual[qual == FALSE] <- "PASS"
  return(qual)
}

filter_df <- marey %>% filter(SEX == "MALE") %>%
    group_by(LG) %>%
    arrange(MB) %>%
    group_by(LG) %>%
    mutate(QC = lgFilter(MB, CM))

filter_df %>%
  ggplot(aes(x = MB, y = CM)) +
  geom_point(aes(color = QC)) +
  facet_wrap(~LG, scales = "free") +
  labs(
    title = "Sex-Averaged markers removed by filtering"
  )

ggsave("sexavg.filtering.png", width = 10, height = 10, units = "in")

filter_df %>% filter(QC == "PASS") %>%
  ggplot(aes(x = MB, y = CM)) +
  geom_point(color = "dodgerblue", shape = 20, alpha = 0.5) +
  facet_wrap(~LG, scales = "free") +
  labs(
    title = "Filtered Sex-Averaged Linkage Maps"
  )
ggsave("sexavg.filtered.png", width = 10, height = 10, units = "in")

marey %>% 
  group_by(LG, SEX) %>%
  arrange(MB) %>%
  group_by(LG) %>%
  mutate(QC = lgFilter(MB, CM)) %>%
  ggplot(aes(x = MB, y = CM)) +
  geom_point(aes(color = QC, shape = SEX)) +
  facet_wrap(~LG, scales = "free") +
  labs(
    title = "Markers removed by filtering"
  )
ggsave("malefemale.filtered.png", width = 10, height = 10, units = "in")
