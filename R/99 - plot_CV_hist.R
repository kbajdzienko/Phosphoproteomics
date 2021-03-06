# Example of usage
# plot_list_hist_cv(file = "CV.pdf", df, df_NA, df_NA_score24, df_NA_score5, df_score24)

#Calculate coeficient of variation within replicates of each group
peakCV <- function(df) {
  df_cv <-
    df$intData %>%
    mutate(intensity = zero.to.na(intensity)) %>%
    left_join(df$sampleData) %>%
    group_by(peak_ID, group) %>%
    summarize(CV = sd(intensity, na.rm = T)/mean(intensity, na.rm = T)) %>%
    spread(group, CV) %>%
    ungroup()
  return(df_cv)
}

#Sum intensities of duplicate phos site entries in each sample
annCV <- function(df) {
  df_cv_sum2 <-
    df$intData %>%
    left_join(select(df$peakData, peak_ID, ann_ID)) %>%
    group_by(ann_ID, sample_ID) %>%
    summarize_at(vars(intensity), sum, na.rm = T) %>%
    mutate(intensity = zero.to.na(intensity)) %>%
    left_join(df$sampleData) %>%
    group_by(ann_ID, group) %>%
    summarize(CV = sd(intensity, na.rm = T)/mean(intensity, na.rm = T)) %>%
    spread(group, CV) %>%
    ungroup()
  return(df_cv_sum2)
}

# Get table with median CV of peak intensities for each sample_group in the experiment
peak_CV_median <- function(df) {
  median.peak.cv <-
    peakCV(df) %>%
    summarize_at(vars(-peak_ID), median, na.rm = T)
  return(median.peak.cv)
}

# Function to plot histograms of several data frames' columns to pdf file
plot_list_hist <- function(df_list, file) {

  # Get the number of columns in dataframes
  n <- unique(sapply(df_list, ncol))
  if (length(n) != 1) stop("Number of columns in data frames is different")

  # Check if all the names are the same between data frames
  name_test <-
    sapply(df_list, function(df) names(df)[-1]) %>%
    apply(1, function(group_ID) length(unique(group_ID)) == 1) %>%
    all()
  if (!name_test) stop("Group names are different in data frames")


  # Open pdf to plot
  pdf(file)
  par(mfrow=c(3,2),
      oma = c(0, 0, 2, 0))

  for (i in 2:n) {
    for (j in 1:length(df_list)) {
      CV <- df_list[[j]][[i]]
      hist(CV,
           xlim = c(0, 2.5),
           breaks = seq(0, 2.5, length.out = 19),
           xlab = "CV",
           main = names(df_list)[j])
      abline(v = median(CV, na.rm = T), col = "red")
    }
    mtext(names(df_list[[1]][i]), outer = TRUE)
    par(mfrow=c(3,2))
  }
  dev.off()
}


# Function to plot CV histograms to pdf file
plot_list_hist_cv <- function(..., file = "OUTPUT/CV_report.pdf") {
  df_names <- as.character(substitute(...()))
  df_list <- c(sapply(list(...), function(df) list(peakCV(df), annCV(df))))
  names(df_list) <- c(sapply(df_names,
                             paste0, c("_cv", "_cv_sum")))
  plot_list_hist(df_list, file)
}

