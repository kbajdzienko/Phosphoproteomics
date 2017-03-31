remove_samples <- function(df, samples){
  #"Removes intensity data of selected samples: remove_samples(df, "23|\\<1\\>|42")
  
  shitsamples <- filter(df$sampleData, grepl(samples, sample_ID)) %>% select(sample_ID)
  df1 <- df
  df1$sampleData <- anti_join(df1$sampleData, shitsamples)
  df1$intData <- anti_join(df1$intData, shitsamples)
  rm(shitsamples)
  return(df1)
}