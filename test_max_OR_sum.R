#Sum or max intensities of duplicate sequence-modification pairs
seqCV <- function(df, method = "sum") {
  df_seq_cv <-
    df$intData %>%
    left_join(select(df$peakData, peak_ID, Sequence, Modifications)) %>%
    unite(SeqMod, Sequence, Modifications) %>%
    group_by(SeqMod, sample_ID) %>%
    summarize_at(vars(intensity), method, na.rm = T) %>%
    mutate(intensity = zero.to.na(intensity)) %>%
    left_join(df$sampleData) %>%
    group_by(SeqMod, group) %>%
    summarize(CV = sd(intensity, na.rm = T)/mean(intensity, na.rm = T)) %>%
    spread(group, CV) %>%
    ungroup()
  return(df_seq_cv)
}

# Median of cv in case of using sum method
cv.sum <- seqCV(df, method = "sum") %>%
  summarize_at(vars(-SeqMod), median, na.rm = T) %>%
  unlist()

# Median of cv in case of using max method
cv.max <- seqCV(df, method = "max") %>%
  summarize_at(vars(-SeqMod), median, na.rm = T) %>%
  unlist()

# Paired statistical test to check whether it's somehow different
wilcox.test(cv.max, cv.sum, paired = T)

# Function to plot CV histograms to pdf file
plot_seq_sum_max_cv <- function(..., file = "max_sum_CV_report.pdf") {
  df_names <- as.character(substitute(...()))
  df_list <- c(sapply(list(...), function(df) list(seqCV(df, "max"), seqCV(df, "sum"))))
  names(df_list) <- c(sapply(df_names,
                             paste0, c("_cv_max", "_cv_sum")))
  plot_list_hist(df_list, file)
}

plot_seq_sum_max_cv(df)
