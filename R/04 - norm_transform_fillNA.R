# Normalization by median intensity of all the peaks in sample
normMedian <- function(df) {
  total.median <- median(df$intData$intensity, na.rm = T)
  df$intData <-
    df$intData %>%
    group_by(sample_ID) %>%
    mutate(intensity = total.median*intensity/median(intensity, na.rm = T)) %>%
    ungroup()
  return(df)
}

# Normalization by mean intensity of all the peaks in sample
normMean <- function(df) {
  total.mean <- mean(df$intData$intensity, na.rm = T)
  df$intData <-
    df$intData %>%
    group_by(sample_ID) %>%
    mutate(intensity = total.mean*intensity/mean(intensity, na.rm = T)) %>%
    ungroup()
  return(df)
}


# Mean-centering and scaling
normScale <- function(df, method = "auto") {
  match.arg(method, c("auto", "pareto"))
  auto <- function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  pareto <- function(x) (x - mean(x, na.rm = T))/sqrt(sd(x, na.rm = T))

  df$intData <-
    df$intData %>%
    group_by(peak_ID) %>%
    mutate_at(vars(intensity), method) %>%
    ungroup()
  return(df)
}


# Fill missing values with PCA-based methods (PPCA or BPCA)
fillNA <- function(df, method = "ppca") {
  match.arg(method, c("ppca", "bpca"))
  if (!any(is.na(df$intData$intensity))) {
    warning("No NA detected")
    return(df)
  }
  intMatrix <-
    intMatrix(df) %>%
    t() %>%
    log10()

  intMatrix[is.na(intMatrix)] <-
    pcaMethods::completeObs(
      pcaMethods::pca(intMatrix, center = T, method = method))[is.na(intMatrix)]

  df$intData <-
    10^intMatrix %>%
    t() %>%
    intData()

  return(df)
}

# Log-transform data
logTransform <- function(df, base = 2) {
  df$intData <- mutate(df$intData, intensity = log(intensity, base = base))
  return(df)
}
