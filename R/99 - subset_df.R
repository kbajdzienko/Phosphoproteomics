subset_df_samples <- function(df, pattern) {
  shitsamples <- filter(df$sampleData, !grepl(pattern, group)) %>% select(sample_ID)
  df$sampleData <- anti_join(df$sampleData, shitsamples)
  df$intData <- anti_join(df$intData, shitsamples)
  df$annIntData <- anti_join(df$annIntData, shitsamples)
  return(df)
}


subset_df_acc <- function(df, acc_list) {
  df$annData <- filter(df$annData, Accession %in% acc_list)
  df$peakData <- filter(df$peakData, Accession %in% acc_list)
  df$intData <- right_join(df$peakData,df$intData, by = "peak_ID") %>% select(peak_ID, sample_ID, intensity)
  df$annIntData <- right_join(df$annData,df$annIntData, by = "ann_ID") %>% select(ann_ID, sample_ID, intensity)
  return(df)
}
