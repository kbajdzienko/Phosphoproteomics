# Filter peaks with too many NA.
# Keep a peak if in at least one group there is less NA, then the set threshold

filter_NA <- function(df, threshold = 0.5) {

  shitpeaks <-
    df$intData %>%
    left_join(df$sampleData) %>%
    mutate(intensity = zero.to.na(intensity)) %>%
    group_by(peak_ID, group) %>%
    summarize(NA_ratio = sum(is.na(intensity))/n()) %>%
    group_by(peak_ID) %>%
    filter(!any(NA_ratio < threshold)) %>%
    ungroup() %>%
    select(peak_ID) %>%
    distinct()

  df$peakData <- anti_join(df$peakData, shitpeaks)
  df$intData <- anti_join(df$intData, shitpeaks)
  return(df)
}


# Filter peaks with too many NA.
# Keep a peak if in at least 2 non-NAs in each group or
# at least 3 non-NAs in 3/4 of groups

filter_NA_Mann <- function(df) {
  shitpeaks <-
    df$intData %>%
    left_join(df$sampleData) %>%
    mutate(intensity = zero.to.na(intensity)) %>%
    group_by(peak_ID, group) %>%
    summarize(non_NA = sum(!is.na(intensity))) %>%
    group_by(peak_ID) %>%
    filter(all(non_NA >= 2) | sum(non_NA >= 3)/n() > 0.75) %>%
    ungroup() %>%
    select(peak_ID) %>%
    distinct()
  df$peakData <- anti_join(df$peakData, shitpeaks)
  df$intData <- anti_join(df$intData, shitpeaks)
  return(df)
}

# Filter peaks with Score less than the set threshold

filter_score <- function(df, score_threshold = 5) {

  shitpeaks <-
    df$peakData %>%
    filter(Score < score_threshold) %>%
    select(peak_ID)

  df$peakData <- anti_join(df$peakData, shitpeaks)
  df$intData <- anti_join(df$intData, shitpeaks)
  return(df)
}

# Filter peaks with MD score less than the set threshold
#MD score refers to confidence of phosphorylation site arrangment

filter_phos_conf <- function(df, score_threshold = 5) {

  shitpeaks <-
    df$peakData %>%
    filter(pep_var_mod_conf < score_threshold) %>%
    select(peak_ID)

  df$peakData <- anti_join(df$peakData, shitpeaks)
  df$intData <- anti_join(df$intData, shitpeaks)
  return(df)
}
