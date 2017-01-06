# Normalization by median intensity of all the peaks in sample
normMedian <- function(df) {
  df$intData <-
    df$intData %>%
    group_by(sample_ID) %>%
    mutate(intensity = intensity/median(intensity, na.rm = T)) %>%
    ungroup()
  return(df)
}

# Fill missing values with PCA-based methods (PPCA or BPCA)
fillNA <- function(df, method = c("ppca", "bpca")) {
  match.arg(method, c("ppca", "bpca"))
  intMatrix <-
    intMatrix(df) %>%
    t() %>%
    log10()

  intMatrix[is.na(intMatrix)] <-
    pcaMethods::completeObs(
      pcaMethods::pca(intMatrix, center = T, method = "ppca"))[is.na(intMatrix)]

  df$intData <-
    10^intMatrix %>%
    t() %>%
    intData()

  return(df)
}

# Fill missing values with PCA-based methods (PPCA or BPCA)
logTransform <- function(df, base = 2) {
  df$intData <- mutate(df$intData, intensity = log(intensity, base = base))
  return(df)
}
