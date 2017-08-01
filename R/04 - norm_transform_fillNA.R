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
normMedian_ann <- function(df) {
  total.median <- median(df$annIntData$intensity, na.rm = T)
  df$annIntData <-
    df$annIntData %>%
    group_by(sample_ID) %>%
    mutate(intensity = total.median*intensity/median(intensity, na.rm = T)) %>%
    ungroup()
  return(df)
}

normMedian2 <- function(df) {
  group_by(df$intData, sample_ID) %>%
    mutate(total = median(intensity, na.rm=T))%>%
    ungroup()%>% mutate(intensity = log2(intensity,na.rm=T), total = log2(total)) %>%
    group_by(sample_ID, peak_ID) %>%
    mutate(intensity = int - total) %>%
    ungroup() %>% mutate(int = 2^int, total=2^total, norm_int=2^norm_int)
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
normScale_ann <- function(df, method = "auto") {
  match.arg(method, c("auto", "pareto"))
  auto <- function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  pareto <- function(x) (x - mean(x, na.rm = T))/sqrt(sd(x, na.rm = T))
  
  df$annIntData <-
    df$annIntData %>%
    group_by(ann_ID) %>%
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

  intMatrix <- t(intMatrix(df))
  if (!is.log(intMatrix)) intMatrix <- log10(intMatrix)

  intMatrix[is.na(intMatrix)] <-
    pcaMethods::completeObs(
      pcaMethods::pca(intMatrix, center = T, method = method))[is.na(intMatrix)]

  if (!is.log(intMatrix(df))) intMatrix <- 10^intMatrix

  df$intData <- intData(t(intMatrix))

  return(df)
}

# Log-transform data
logTransform <- function(df, base = 2) {
  if (any(df$intData$intensity < 0, na.rm = T)) {
    warning("Negative values detected.\nData are probably already log-transformed or scaled.\nNo transformation applied.")
    return(df)
  }
  if (is.log(df$intData$intensity)) {
    warning("Data are already log-transformed, no transformation applied.")
    return(df)
  }
  df$intData <- mutate(df$intData, intensity = log(intensity, base = base))
  if (!is.null(df$annData))
    df$annIntData <- mutate(df$annIntData, intensity = log(intensity, base = base))
  return(df)
}

logTransform_ann <- function(df, base = 2) {
  if (any(df$annIntData$intensity < 0, na.rm = T)) {
    warning("Negative values detected.\nData are probably already log-transformed or scaled.\nNo transformation applied.")
    return(df)
  }
  if (is.log(df$annIntData$intensity)) {
    warning("Data are already log-transformed, no transformation applied.")
    return(df)
  }
  
    df$annIntData <- mutate(df$annIntData, intensity = log(intensity, base = base))
  return(df)
}
