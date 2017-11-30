remove_samples <- function(df, sample_pattern){
  #"Removes intensity data of selected samples: remove_samples(df, "23|\\<1\\>|42")
  
  shitsamples <- filter(df$sampleData, grepl(sample_pattern, sample_ID)) %>% select(sample_ID)
  df1 <- df
  df1$sampleData <- anti_join(df1$sampleData, shitsamples)
  df1$intData <- anti_join(df1$intData, shitsamples)
  rm(shitsamples)
  return(df1)
}

remove_samples2 <- function(df, sample_vector){
  #"Removes intensity data of selected samples: remove_samples(df, c("1","2"))
  
  shitsamples <- filter(df$sampleData, sample_ID %in% sample_vector) %>% select(sample_ID)
  df1 <- df
  df1$sampleData <- anti_join(df1$sampleData, shitsamples)
  df1$intData <- anti_join(df1$intData, shitsamples)
  df1$annIntData <- anti_join(df1$annIntData, shitsamples)
  rm(shitsamples)
  return(df1)
}
