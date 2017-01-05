# Sum intensities for the same protein, the same modification
sumANN <- function(df) {
  df$intData <-
    df$intData %>%
    left_join(select(df$peakData, peak_ID, ann_ID)) %>%
    group_by(ann_ID, sample_ID) %>%
    summarize_at(vars(intensity), sum, na.rm = T) %>%
    mutate(intensity = zero.to.na(intensity)) %>%
    ungroup()
  return(df)
}

# Normalization by median intensity of all the peaks in sample
normMedian <- function(df) {
  df$intData <-
    df$intData %>%
    group_by(sample_ID) %>%
    mutate(intensity = intensity/median(intensity, na.rm = T)) %>%
    ungroup()
  return(df)
}

# Fill missing values
fillNA <- function(df, method = c("ppca", "bpca")) {
  match.arg(method, c("ppca", "bpca"))
  intMatrix <- intTable(df) %>%
    as.data.frame() %>%
    'row.names<-'(.$peak_ID) %>% head()
    select(-peak_ID) %>%
    t() %>%
    log10()

  intMatrix[is.na(intMatrix)] <-
    pcaMethods::completeObs(
      pcaMethods::pca(intMatrix, center = T, method = "ppca"))[is.na(intMatrix)]

  df$intData <-
    10^intMatrix %>%
    t() %>%
    as.data.frame() %>%
    mutate(peak_ID = row.names(.)) %>%
    select(peak_ID, everything()) %>%
    gather("sample_ID", "intensity", -peak_ID) %>%
    tbl_df()

  return(df)
}

logTransform <- function(df) {
  df$intData <- mutate_at(df$intData, vars(intensity), log10)
  return(df)
}
