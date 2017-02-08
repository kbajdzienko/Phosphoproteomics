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
  df$intData <- semi_join(df$peakData, by = Accession)
  df$annIntData <- semi_join(df$annData, by = Accession)
  return(df)
}
